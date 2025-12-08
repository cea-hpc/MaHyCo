// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#include "PerfectGasEOSService.h"
#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"

#include <accenv/IAccEnv.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Constructeur de la classe                                                 */
/*---------------------------------------------------------------------------*/
PerfectGasEOSService::PerfectGasEOSService(const ServiceBuildInfo & sbi)
  : ArcanePerfectGasEOSObject(sbi) {
  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::initEOS(IMeshEnvironment* env)
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
/* Formule PerfectGas pour calculer unitairement pression, vitesse du son et dp/de */
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

void PerfectGasEOSService::applyEOS(IMeshEnvironment* env)
{
  PROF_ACC_BEGIN(__FUNCTION__);
  Real adiabatic_cst = getAdiabaticCst(env);
  // Calcul de la pression et de la vitesse du son
#if 0
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;
    if (m_density[ev] == 0.) info() << ev.globalCell().localId() << " densité nulle";
    compute_pressure_sndspd_PG(adiabatic_cst,
        m_density[ev], m_internal_energy[ev],
        m_pressure[ev], m_sound_speed[ev], m_dpde[ev]);
  }
#elif 0
  Parallel::Foreach(env->envView(),[&](EnvItemVectorView view)
  {
    ENUMERATE_ENVCELL(ienvcell,view){
      compute_pressure_sndspd_PG(adiabatic_cst,
          m_density[ienvcell], m_internal_energy[ienvcell],
          m_pressure[ienvcell], m_sound_speed[ienvcell], m_dpde[ienvcell]);
    }
  });
#else

  auto queue = m_acc_env->newQueue();
  {
    auto command = makeCommand(queue);

    auto in_density         = ax::viewIn (command, m_density);
    auto in_internal_energy = ax::viewIn (command, m_internal_energy);

    auto out_pressure       = ax::viewOut(command, m_pressure);
    auto out_sound_speed    = ax::viewOut(command, m_sound_speed);
    auto out_dpde           = ax::viewOut(command, m_dpde);

    command << RUNCOMMAND_MAT_ENUMERATE(EnvCell, evi, env) {

      Real pressure, sound_speed, dpde;

      compute_pressure_sndspd_PG(adiabatic_cst,
          in_density[evi], in_internal_energy[evi],
          pressure, sound_speed, dpde);

      out_pressure[evi] = pressure;
      out_sound_speed[evi] = sound_speed;
      out_dpde[evi] = dpde;

    };
  }

#endif
  PROF_ACC_END;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  Real adiabatic_cst = getAdiabaticCst(env);
  // Calcul de la pression et de la vitesse du son
#if 0
    Real internal_energy = m_internal_energy[ev];
    Real density = m_density[ev];
    if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
    Real pressure = (adiabatic_cst - 1.) * density * internal_energy;
    m_pressure[ev] = pressure;
    m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
    m_dpde[ev] = (adiabatic_cst - 1.) * density;
#else
    if (m_density[ev] == 0.) info() << ev.globalCell().localId() << " densité nulle";
    compute_pressure_sndspd_PG(adiabatic_cst,
        m_density[ev], m_internal_energy[ev],
        m_pressure[ev], m_sound_speed[ev], m_dpde[ev]);
#endif
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real PerfectGasEOSService::getAdiabaticCst([[maybe_unused]] IMeshEnvironment* env) { return options()->adiabaticCst();}
Real PerfectGasEOSService::getTensionLimitCst([[maybe_unused]] IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_PERFECTGASEOS(PerfectGas, PerfectGasEOSService);
