// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "StdPerfectGasAcc1EOSService.h"
#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemVector.h"

#include <accenv/IAccEnv.h>
#include <eos/stdperfectgasacc1/EosTypes.h>
#include <eos/stdperfectgasacc1/PhyVarType.h>

using namespace Arcane;
using namespace Arcane::Materials;

namespace Stdperfectgasacc1 {

/*---------------------------------------------------------------------------*/
/* Constructeur de la classe                                                 */
/*---------------------------------------------------------------------------*/
StdPerfectGasAcc1EOSService::StdPerfectGasAcc1EOSService(const ServiceBuildInfo & sbi)
  : ArcaneStdPerfectGasAcc1EOSObject(sbi) {
  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasAcc1EOSService::initEOS(IMeshEnvironment* env)
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
/* Formule StdPerfectGasAcc1 pour calculer unitairement pression, vitesse du son et dp/de */
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

void add_phy_var(std::vector<PhyVarType>& phy_vars,
    Arcane::Materials::MaterialVariableCellReal var, 
    const Arcane::Materials::MatCellVectorView* matcell_vector_view) 
{
  prof_acc_begin(var.name().localstr());
  phy_vars.emplace_back(var, matcell_vector_view);
  prof_acc_end(var.name().localstr());
}

void phy_var_2_arc_var(const PhyVarType& phy_var,
    Arcane::Materials::MaterialVariableCellReal var,
    const Arcane::Materials::MatCellVectorView* matcell_vector_view)
{
  prof_acc_begin(var.name().localstr());
  phy_var.copyRawDataToVar(var, matcell_vector_view);
  prof_acc_end(var.name().localstr());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasAcc1EOSService::applyEOS(IMeshEnvironment* env)
{
  prof_acc_begin(__FUNCTION__);
  // Calcul de la pression et de la vitesse du son

  ENUMERATE_MAT(mat_i, env) 
  {
    IMeshMaterial* mat = *mat_i;

    // std_matcell_vector construction
    prof_acc_begin("matcell_vector_view"); 
    Arcane::Materials::MatCellVectorView matcell_vector_view(mat->matView());
    prof_acc_end("matcell_vector_view");

    // Arcane var => PhyVarTpe
    prof_acc_begin("ArcVar=>PhyVarTpe"); 
    std::vector<PhyVarType> phy_vars;

    add_phy_var(phy_vars, m_density        , &matcell_vector_view);
    add_phy_var(phy_vars, m_internal_energy, &matcell_vector_view);
    add_phy_var(phy_vars, m_pressure       , &matcell_vector_view);
    add_phy_var(phy_vars, m_sound_speed    , &matcell_vector_view);
    add_phy_var(phy_vars, m_dpde           , &matcell_vector_view);
    prof_acc_end("ArcVar=>PhyVarTpe");

    // The EOS itself
    prof_acc_begin("EOS_PG"); 
    Real adiabatic_cst = getAdiabaticCst(env);

    auto queue = m_acc_env->newQueue(); // synchronous queue
    auto command = makeCommand(queue);

    auto nelt = phy_vars[0].raw_data.size();

    Span<const Real> in_density         (phy_vars[0].raw_data.data(), nelt);
    Span<const Real> in_internal_energy (phy_vars[1].raw_data.data(), nelt);

    Span<Real>       out_pressure       (phy_vars[2].raw_data.data(), nelt);
    Span<Real>       out_sound_speed    (phy_vars[3].raw_data.data(), nelt);
    Span<Real>       out_dpde           (phy_vars[4].raw_data.data(), nelt);

    command << RUNCOMMAND_LOOP1(iter, nelt)
    {
      auto i = iter()[0];

      compute_pressure_sndspd_PG(adiabatic_cst,
	  in_density[i], in_internal_energy[i],
	  out_pressure[i], out_sound_speed[i], out_dpde[i]);
    };
    prof_acc_end("EOS_PG");

    // PhyVarType => Arcane var
    prof_acc_begin("PhyVarTpe=>ArcVar"); 
    phy_var_2_arc_var(phy_vars[2], m_pressure    , &matcell_vector_view);
    phy_var_2_arc_var(phy_vars[3], m_sound_speed , &matcell_vector_view);
    phy_var_2_arc_var(phy_vars[4], m_dpde        , &matcell_vector_view);
    prof_acc_end("PhyVarTpe=>ArcVar");
  }

  prof_acc_end(__FUNCTION__);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasAcc1EOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
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
Real StdPerfectGasAcc1EOSService::getAdiabaticCst([[maybe_unused]] IMeshEnvironment* env) { return options()->adiabaticCst();}
Real StdPerfectGasAcc1EOSService::getTensionLimitCst([[maybe_unused]] IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_STDPERFECTGASACC1EOS(StdPerfectGasAcc1, StdPerfectGasAcc1EOSService);

} // Stdperfectgasacc1

