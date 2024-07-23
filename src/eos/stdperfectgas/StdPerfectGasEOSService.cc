// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#include "StdPerfectGasEOSService.h"
#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemVector.h"

#include <accenv/IAccEnv.h>
#include <eos/stdperfectgas/EosTypes.h>
#include <eos/stdperfectgas/PhyVarType.h>

using namespace Arcane;
using namespace Arcane::Materials;

namespace Stdperfectgas {

/*---------------------------------------------------------------------------*/
/* Constructeur de la classe                                                 */
/*---------------------------------------------------------------------------*/
StdPerfectGasEOSService::StdPerfectGasEOSService(const ServiceBuildInfo & sbi)
  : ArcaneStdPerfectGasEOSObject(sbi) {
  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasEOSService::initEOS(IMeshEnvironment* env)
{
  Real adiabatic_cst = getAdiabaticCst(env);
  // Initialise l'énergie et la vitesse du son
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real pressure = m_pressure[ev];
    Real density = m_density[ev];
    m_internal_energy[ev] = pressure / ((adiabatic_cst - 1.) * density);
    m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
  }
}

/*---------------------------------------------------------------------------*/
/* Formule StdPerfectGas pour calculer unitairement pression, vitesse du son et dp/de */
/* Cette formule est appelée dans applyEOS(...) et applyOneCellEOS(...)      */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE inline void compute_pressure_sndspd_PG(Real adiabatic_cst,
    Real density, Real internal_energy,
    Real& pressure, Real& sound_speed, Real& dpde) 
{
  pressure = (adiabatic_cst - 1.) * density * internal_energy;
  sound_speed = sqrt(adiabatic_cst * pressure / density);
  dpde = (adiabatic_cst - 1.) * density;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasEOSService::applyEOS(IMeshEnvironment* env)
{
  prof_acc_begin(__FUNCTION__);
  // Calcul de la pression et de la vitesse du son

  prof_acc_begin("cell_vector"); 
  CellVector cell_vector(mesh()->cellFamily());
  ENUMERATE_ENVCELL(ienvcell,env->envView())
  {
    cell_vector.addItem((*ienvcell).globalCell());
  }
  prof_acc_end("cell_vector");

  ENUMERATE_MAT(mat_i, env) 
  {
    IMeshMaterial* mat = *mat_i;

    // std_matcell_vector construction
    prof_acc_begin("std_matcell_vector"); 
    StdMatCellVector std_matcell_vector;
    std_matcell_vector.reserve(mat->matView().nbItem());
    MatCellVector mat_cells(cell_vector.view(), mat);
    ENUMERATE_MATCELL(matcell_i, mat_cells.view())
    {
      std_matcell_vector.push_back(*matcell_i);
    }
    prof_acc_end("std_matcell_vector");

    // Arcane var => PhyVarTpe
    prof_acc_begin("ArcVar=>PhyVarTpe"); 
    std::vector<PhyVarType> phy_vars;

    phy_vars.emplace_back(m_density        , &std_matcell_vector);
    phy_vars.emplace_back(m_internal_energy, &std_matcell_vector);
    phy_vars.emplace_back(m_pressure       , &std_matcell_vector);
    phy_vars.emplace_back(m_sound_speed    , &std_matcell_vector);
    phy_vars.emplace_back(m_dpde           , &std_matcell_vector);
    prof_acc_end("ArcVar=>PhyVarTpe");

    // The EOS itself
    prof_acc_begin("EOS_PG"); 
    Real adiabatic_cst = getAdiabaticCst(env);

    auto& density         (phy_vars[0].raw_data);
    auto& internal_energy (phy_vars[1].raw_data);
    auto& pressure        (phy_vars[2].raw_data);
    auto& sound_speed     (phy_vars[3].raw_data);
    auto& dpde            (phy_vars[4].raw_data);

    for(size_t i(0) ; i < phy_vars[0].raw_data.size() ; ++i) 
    {
      compute_pressure_sndspd_PG(adiabatic_cst,
	  density[i], internal_energy[i],
	  pressure[i], sound_speed[i], dpde[i]);
    }
    prof_acc_end("EOS_PG");

    // PhyVarType => Arcane var
    prof_acc_begin("PhyVarTpe=>ArcVar"); 
    phy_vars[2].copyRawDataToVar(m_pressure    , &std_matcell_vector);
    phy_vars[3].copyRawDataToVar(m_sound_speed , &std_matcell_vector);
    phy_vars[4].copyRawDataToVar(m_dpde        , &std_matcell_vector);
    prof_acc_end("PhyVarTpe=>ArcVar");
  }

  prof_acc_end(__FUNCTION__);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  Real adiabatic_cst = getAdiabaticCst(env);
  // Calcul de la pression et de la vitesse du son
  if (m_density[ev] == 0.) info() << ev.globalCell().localId() << " densité nulle";
  compute_pressure_sndspd_PG(adiabatic_cst,
      m_density[ev], m_internal_energy[ev],
      m_pressure[ev], m_sound_speed[ev], m_dpde[ev]);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real StdPerfectGasEOSService::getAdiabaticCst([[maybe_unused]] IMeshEnvironment* env) { return options()->adiabaticCst();}
Real StdPerfectGasEOSService::getTensionLimitCst([[maybe_unused]] IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_STDPERFECTGASEOS(StdPerfectGas, StdPerfectGasEOSService);

} // Stdperfectgas

