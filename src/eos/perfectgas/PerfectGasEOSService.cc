// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "PerfectGasEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::initEOS(IMeshEnvironment* env)
{
  Real adiabatic_cst = options()->adiabaticCst();
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
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::applyEOS(IMeshEnvironment* env)
{
  Real adiabatic_cst = options()->adiabaticCst();
  // Calcul de la pression et de la vitesse du son
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real internal_energy = m_internal_energy[ev];
    Real density = m_density[ev];
    Real pressure = (adiabatic_cst - 1.) * density * internal_energy;
    m_pressure[ev] = pressure;
    Cell cell = ev.globalCell();
    if (pressure < 0.) info() << cell.localId() << " " << pressure << " " << internal_energy << " " << density << " et adia " << adiabatic_cst;
    m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real PerfectGasEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_PERFECTGASEOS(PerfectGas, PerfectGasEOSService);
