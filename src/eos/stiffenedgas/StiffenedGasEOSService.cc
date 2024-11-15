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
     // Cv et indice adiabatique
  Real cv  = options()->specificHeat;
  Real gamma = options()->adiabaticCst;
  Real limit_tension = options()->limitTension;
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real pressure = m_pressure[ev];
    Real density = m_density[ev];
    m_internal_energy[ev] = (pressure + (gamma * limit_tension)) / ((gamma - 1.) * density);
    m_sound_speed[ev] = sqrt((gamma/density)*(pressure+limit_tension));
    m_density_0[ev] = m_density[ev];
    // calcul de la temperature de la constante (on prend celle de l'air : 287.)
    m_temperature[ev] = pressure / (287. * density);
    m_internal_energy_0[ev] =  m_internal_energy[ev];
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::ReinitEOS(IMeshEnvironment* env)
{
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::applyEOS(IMeshEnvironment* env)
{
  // Cv et indice adiabatique
  Real cv  = options()->specificHeat;
  Real gamma = options()->adiabaticCst;
  Real limit_tension = options()->limitTension;
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    if (m_maille_endo[ev.globalCell()] == 0) {
        Real internal_energy = m_internal_energy[ev];
        Real density = m_density[ev];
        Real pressure = ((gamma - 1.) * density * internal_energy) - (gamma * limit_tension);
        m_pressure[ev] = pressure;
        m_sound_speed[ev] = sqrt((gamma/density)*(pressure+limit_tension));
        m_dpde[ev] = (gamma - 1.) * density;
        m_temperature_n[ev] = m_temperature[ev];
        m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / cv + m_temperature_n[ev];
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  if (m_maille_endo[ev.globalCell()] == 1) return;
  
  // Cv et indice adiabatique
  Real cv  = options()->specificHeat;
  Real gamma = options()->adiabaticCst;
  Real limit_tension = options()->limitTension;
  Real internal_energy = m_internal_energy[ev];
  Real density = m_density[ev];
  Real pressure = ((gamma - 1.) * density * internal_energy) - (gamma * limit_tension);
  m_pressure[ev] = pressure;
  m_sound_speed[ev] = sqrt((gamma/density)*(pressure+limit_tension));
  m_dpde[ev] = (gamma - 1.) * density;
        m_temperature_n[ev] = m_temperature[ev];
  m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / cv + m_temperature_n[ev];
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::Endommagement(IMeshEnvironment* env)
{
  Real damage_thresold = options()->tensionDamageThresold();
  Real density_thresold = options()->densityDamageThresold();
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;  
    Cell cell = ev.globalCell();  
    if (m_maille_endo[ev.globalCell()] == 0) {
        // Maille saine : verification des seuils 
        if (m_pressure[cell] < damage_thresold || m_density[cell]/m_density_0[cell] < density_thresold) {
            // maille devient endommagée
            m_maille_endo[ev] = 1;
            m_density_fracture[ev] = m_density[ev];
            m_internal_energy_fracture[ev] = m_internal_energy[ev];
            m_pressure[ev] = 0.;
        }
    } else { 
        // Maille endo : on verifie si elle ne s'est pas recompactée
        if (m_density[ev] > m_density_fracture[ev]) {
            // c'est le cas : on la déclare saine
            m_maille_endo[ev] = 0;
            m_internal_energy[ev] = m_internal_energy_fracture[ev];
        }
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real StiffenedGasEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real StiffenedGasEOSService::getdensityDamageThresold(IMeshEnvironment* env) { return options()->densityDamageThresold();}
/*--------------------------------------------- */
/*--------------------------------------------- */
ARCANE_REGISTER_SERVICE_STIFFENEDGASEOS(StiffenedGas, StiffenedGasEOSService);
