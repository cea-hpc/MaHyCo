﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "StiffenedGasEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::initEOS(IMeshEnvironment* env)
{
  // Initialise l'énergie et la vitesse du son
  Real limit_tension = options()->limitTension();
  Real adiabatic_cst = options()->adiabaticCst();
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
  Real limit_tension = options()->limitTension();
  Real adiabatic_cst = options()->adiabaticCst();
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real internal_energy = m_internal_energy[ev];
    Real density = m_density[ev];
    Real pressure = ((adiabatic_cst - 1.) * density * internal_energy) - (adiabatic_cst * limit_tension);
    m_pressure[ev] = pressure;
    m_sound_speed[ev] = sqrt((adiabatic_cst/density)*(pressure+limit_tension));
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real StiffenedGasEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_STIFFENEDGASEOS(StiffenedGas, StiffenedGasEOSService);
