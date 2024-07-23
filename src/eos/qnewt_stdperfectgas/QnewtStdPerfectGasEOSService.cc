// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#include "QnewtStdPerfectGasEOSService.h"
#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemVector.h"

#include <accenv/IAccEnv.h>
#include <eos/qnewt_stdperfectgas/EosTypes.h>
#include <eos/qnewt_stdperfectgas/PhyVarType.h>
#include <eos/qnewt_stdperfectgas/Utils.h>

using namespace Arcane;
using namespace Arcane::Materials;

namespace Qnewt_stdperfectgas {

/*---------------------------------------------------------------------------*/
/* Constructeur de la classe                                                 */
/*---------------------------------------------------------------------------*/
QnewtStdPerfectGasEOSService::QnewtStdPerfectGasEOSService(const ServiceBuildInfo & sbi)
  : ArcaneQnewtStdPerfectGasEOSObject(sbi) {
  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void QnewtStdPerfectGasEOSService::initEOS(IMeshEnvironment* env)
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
/* Formule QnewtStdPerfectGas pour calculer unitairement pression, vitesse du son et dp/de */
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

void QnewtStdPerfectGasEOSService::applyEOS(IMeshEnvironment* env)
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

    prof_acc_begin("Newton");

    prof_acc_begin("not_conv_matcell_vector"); 
    StdMatCellVector not_conv_matcell_vector(std_matcell_vector);
    prof_acc_end("not_conv_matcell_vector"); 

    Integer newt_max = 30;
    bool remain_item = true;
    Integer inewt;
    for(inewt = 0 ; inewt<newt_max && remain_item ; ++inewt)
    {
      prof_acc_begin("IterNewt");
      // Direct call
      {
	// Arcane var => PhyVarTpe
	prof_acc_begin("ArcVar=>PhyVarTpe"); 
	std::vector<PhyVarType> phy_vars;

	phy_vars.emplace_back(m_density        , &not_conv_matcell_vector);
	phy_vars.emplace_back(m_internal_energy, &not_conv_matcell_vector);
	phy_vars.emplace_back(m_pressure       , &not_conv_matcell_vector);
	phy_vars.emplace_back(m_sound_speed    , &not_conv_matcell_vector);
	phy_vars.emplace_back(m_dpde           , &not_conv_matcell_vector);
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
	phy_vars[2].copyRawDataToVar(m_pressure    , &not_conv_matcell_vector);
	phy_vars[3].copyRawDataToVar(m_sound_speed , &not_conv_matcell_vector);
	phy_vars[4].copyRawDataToVar(m_dpde        , &not_conv_matcell_vector);
	prof_acc_end("PhyVarTpe=>ArcVar");
      }

      prof_acc_begin("to_erase");
      Integer nb_elem = not_conv_matcell_vector.size();
#if 0
      Integer nb_var=1;
      std::vector<std::vector<Integer>> conv_item(nb_var, std::vector<Integer>(nb_elem, 0));
#endif
      UniqueArray<bool> to_erase(nb_elem, false);
      for(Integer ielem=0 ; ielem<nb_elem ; ++ielem) {
	to_erase[ielem] = (ielem==0) || (ielem%10);
      }
      prof_acc_end("to_erase");

      prof_acc_begin("resizeFromMask");
      resizeObjectFromMask(to_erase, &(not_conv_matcell_vector));
      prof_acc_end("resizeFromMask");

#if 0
      for(Integer iva=0 ; ivar<nb_var ; ++ivar) {
	resizeObjectFromMask(to_erase, &(conv_item[ivar]));
      }
#endif

      // Stop criteria
      remain_item = !(not_conv_matcell_vector.empty());
      /*
      info() << "iter=" << inewt 
	<< ", to_erase.size()=" << std::count_if(to_erase.begin(), to_erase.end(), [](bool e) {return e;} ) 
	<< ", nremain=" << not_conv_matcell_vector.size();
       */
      prof_acc_end("IterNewt");
    }
    //info() << "Number of Newton iterations=" << inewt+1;
    prof_acc_end("Newton"); 
  }

  prof_acc_end(__FUNCTION__);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void QnewtStdPerfectGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
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
Real QnewtStdPerfectGasEOSService::getAdiabaticCst([[maybe_unused]] IMeshEnvironment* env) { return options()->adiabaticCst();}
Real QnewtStdPerfectGasEOSService::getTensionLimitCst([[maybe_unused]] IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_QNEWTSTDPERFECTGASEOS(QnewtStdPerfectGas, QnewtStdPerfectGasEOSService);

} // Qnewt_stdperfectgas

