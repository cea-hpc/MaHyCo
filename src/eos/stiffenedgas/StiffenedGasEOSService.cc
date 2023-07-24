// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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

void StiffenedGasEOSService::ReinitEOS(IMeshEnvironment* env)
{
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
    if (m_maille_endo[ev] == 0) {
        Real internal_energy = m_internal_energy[ev];
        Real density = m_density[ev];
        Real pressure = ((adiabatic_cst - 1.) * density * internal_energy) - (adiabatic_cst * limit_tension);
        m_pressure[ev] = pressure;
        m_sound_speed[ev] = sqrt((adiabatic_cst/density)*(pressure+limit_tension));
        m_dpde[ev] = (adiabatic_cst - 1.) * density;
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void StiffenedGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  if (m_maille_endo[ev] == 1) return;
  
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

void StiffenedGasEOSService::Endommagement(IMeshEnvironment* env)
{
  Real damage_thresold = options()->damageThresold();
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    if (m_maille_endo[ev] == 0) {
        // Maille saine : verification des seuils 
        if (m_pressure[ev] < damage_thresold) {
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
Real StiffenedGasEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real StiffenedGasEOSService::getSpecificHeatCst(IMeshEnvironment* env) { return options()->specificHeat();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_STIFFENEDGASEOS(StiffenedGas, StiffenedGasEOSService);
