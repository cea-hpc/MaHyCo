// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "PerfectGasEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::initEOS(IMeshEnvironment* env)
{
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  // Initialise l'énergie et la vitesse du son pour chaque maille de l'environnement
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real pressure = m_pressure[ev];
    Real density = m_density[ev];
    Cell cell = ev.globalCell();
    // Affiche l'identifiant local de la maille, la pression et la densité
    m_internal_energy[ev] = pressure / ((adiabatic_cst - 1.) * density);
    m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
    // calcul de la temperature en fonction de la chaleur specifique
    m_temperature[ev] = m_internal_energy[ev] / specific_heat;
    
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::ReinitEOS(IMeshEnvironment* env)
{
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::applyEOS(IMeshEnvironment* env)
{
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  // Calcul de la pression et de la vitesse du son pour chaque maille de l'environnement
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;
    if (m_maille_endo[ev.globalCell()] == 0) { 
        Real internal_energy = m_internal_energy[ev];
        Real density = m_density[ev];
        if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
        Real pressure = (adiabatic_cst - 1.) * density * internal_energy;
        m_pressure[ev] = pressure;
        m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
        m_dpde[ev] = (adiabatic_cst - 1.) * density;
        // calcul de la temperature en fonction de la chaleur specifique
        m_temperature[ev] = internal_energy / specific_heat;
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  if (m_maille_endo[ev.globalCell()] == 1) return;
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  // Calcul de la pression,la vitesse du son  et le gradient de pression pour une maille donnée
  Real internal_energy = m_internal_energy[ev];
  Real density = m_density[ev];
  if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
  Real pressure = (adiabatic_cst - 1.) * density * internal_energy;
  m_pressure[ev] = pressure;
  m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
  m_dpde[ev] = (adiabatic_cst - 1.) * density;
  // calcul de la temperature en fonction de la chaleur specifique
  m_temperature[ev] = internal_energy / specific_heat;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::Endommagement(IMeshEnvironment* env)
{
  Real damage_thresold = options()->damageThresold();
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    if (m_maille_endo[ev.globalCell()] == 0) {
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
Real PerfectGasEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real PerfectGasEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real PerfectGasEOSService::getSpecificHeatCst(IMeshEnvironment* env) { return options()->specificHeat();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_PERFECTGASEOS(PerfectGas, PerfectGasEOSService);
