// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "StdPerfectGasAcc2EOSService.h"
#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemVector.h"

#include <accenv/IAccEnv.h>
#include <eos/stdperfectgasacc2/EosTypes.h>
#include <eos/stdperfectgasacc2/PhyVars.h>

namespace Stdperfectgasacc2 {

/*---------------------------------------------------------------------------*/
/* Constructeur de la classe                                                 */
/*---------------------------------------------------------------------------*/
StdPerfectGasAcc2EOSService::StdPerfectGasAcc2EOSService(const ServiceBuildInfo & sbi)
  : ArcaneStdPerfectGasAcc2EOSObject(sbi) {
  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasAcc2EOSService::initEOS(IMeshEnvironment* env)
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
/* Formule StdPerfectGasAcc2 pour calculer unitairement pression, vitesse du son et dp/de */
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

void StdPerfectGasAcc2EOSService::applyEOS(IMeshEnvironment* env)
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
    m_phy_vars.addPhyVar(m_density        , &matcell_vector_view);
    m_phy_vars.addPhyVar(m_internal_energy, &matcell_vector_view);
    m_phy_vars.addPhyVar(m_pressure       , &matcell_vector_view);
    m_phy_vars.addPhyVar(m_sound_speed    , &matcell_vector_view);
    m_phy_vars.addPhyVar(m_dpde           , &matcell_vector_view);
    prof_acc_end("ArcVar=>PhyVarTpe");

    //
    m_phy_vars.copyFromAllArcVars();

    // The EOS itself
    prof_acc_begin("EOS_PG"); 
    Real adiabatic_cst = getAdiabaticCst(env);

    auto queue = m_acc_env->newQueue(); // synchronous queue
    auto command = makeCommand(queue);

    auto nelt = m_phy_vars.rawData(0).size();

    Span<const Real> in_density         (m_phy_vars.rawData(0));
    Span<const Real> in_internal_energy (m_phy_vars.rawData(1));
    Span<Real>       out_pressure       (m_phy_vars.rawData(2));
    Span<Real>       out_sound_speed    (m_phy_vars.rawData(3));
    Span<Real>       out_dpde           (m_phy_vars.rawData(4));

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
    m_phy_vars.copyIntoArcVar(2);
    m_phy_vars.copyIntoArcVar(3);
    m_phy_vars.copyIntoArcVar(4);
    prof_acc_end("PhyVarTpe=>ArcVar");

    m_phy_vars.clear();
  }

  prof_acc_end(__FUNCTION__);
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StdPerfectGasAcc2EOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
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
Real StdPerfectGasAcc2EOSService::getAdiabaticCst([[maybe_unused]] IMeshEnvironment* env) { return options()->adiabaticCst();}
Real StdPerfectGasAcc2EOSService::getTensionLimitCst([[maybe_unused]] IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_STDPERFECTGASACC2EOS(StdPerfectGasAcc2, StdPerfectGasAcc2EOSService);

} // Stdperfectgasacc2

