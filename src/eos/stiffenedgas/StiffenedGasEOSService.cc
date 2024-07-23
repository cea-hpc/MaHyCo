// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#include "StiffenedGasEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::initEOS(IMeshEnvironment* env)
{
  // Initialise l'énergie et la vitesse du son
  Real limit_tension = getTensionLimitCst(env);
  Real adiabatic_cst = getAdiabaticCst(env);
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real pressure = m_pressure[ev];
    Real density = m_density[ev];
    m_internal_energy[ev] = (pressure + (adiabatic_cst * limit_tension)) / ((adiabatic_cst - 1.) * density);
    m_sound_speed[ev] = sqrt((adiabatic_cst/density)*(pressure+limit_tension));
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::applyEOS(IMeshEnvironment* env)
{
  // Calcul de la pression et de la vitesse du son
  Real limit_tension = getTensionLimitCst(env);
  Real adiabatic_cst = getAdiabaticCst(env);
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real internal_energy = m_internal_energy[ev];
    Real density = m_density[ev];
    Real pressure = ((adiabatic_cst - 1.) * density * internal_energy) - (adiabatic_cst * limit_tension);
    m_pressure[ev] = pressure;
    m_sound_speed[ev] = sqrt((adiabatic_cst/density)*(pressure+limit_tension));
    m_dpde[ev] = (adiabatic_cst - 1.) * density;
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  // Calcul de la pression et de la vitesse du son
  Real limit_tension = getTensionLimitCst(env);
  Real adiabatic_cst = getAdiabaticCst(env);
    Real internal_energy = m_internal_energy[ev];
    Real density = m_density[ev];
    Real pressure = ((adiabatic_cst - 1.) * density * internal_energy) - (adiabatic_cst * limit_tension);
    m_pressure[ev] = pressure;
    m_sound_speed[ev] = sqrt((adiabatic_cst/density)*(pressure+limit_tension));
    m_dpde[ev] = (adiabatic_cst - 1.) * density;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real StiffenedGasEOSService::getAdiabaticCst([[maybe_unused]] IMeshEnvironment* env) { return options()->adiabaticCst();}
Real StiffenedGasEOSService::getTensionLimitCst([[maybe_unused]] IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_STIFFENEDGASEOS(StiffenedGas, StiffenedGasEOSService);
