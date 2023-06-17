// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MieGruneisenEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MieGruneisenEOSService::initEOS(IMeshEnvironment* env)
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
    Real J = options()->rho_0 / density ;
    Real mu = 1./J -1.;
    // Affiche l'identifiant local de la maille, la pression et la densité
    m_internal_energy[ev] = (pressure - options()->C * mu + options()->D * pow(mu,2.) + options()->S * pow(mu,3)) / (adiabatic_cst * density);
    m_sound_speed[ev] = options()->C / options()->rho_0;
    // calcul de la temperature en fonction de la chaleur specifique
    m_temperature[ev] = m_internal_energy[ev] / specific_heat;
    
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MieGruneisenEOSService::applyEOS(IMeshEnvironment* env)
{
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  // Calcul de la pression et de la vitesse du son pour chaque maille de l'environnement
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real internal_energy = m_internal_energy[ev];
    Real density = m_density[ev];
    if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
    Real J = options()->rho_0 / density ;
    Real mu = 1./J -1.;
    m_pressure[ev] = options()->C * mu + options()->D * pow(mu,2.) + options()->S * pow(mu,3) +   adiabatic_cst * density * internal_energy;
    m_temperature[ev] = internal_energy / specific_heat;
    m_sound_speed[ev] = options()->C / options()->rho_0;
    m_dpde[ev] = 0.;
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MieGruneisenEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  // Calcul de la pression,la vitesse du son  et le gradient de pression pour une maille donnée
  Real internal_energy = m_internal_energy[ev];
  Real density = m_density[ev];
  if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
  Real J = options()->rho_0 / density ;
    Real mu = 1./J -1.;
    m_pressure[ev] = options()->C * mu + options()->D * pow(mu,2.) + options()->S * pow(mu,3) +   adiabatic_cst * density * internal_energy;
    m_sound_speed[ev] = options()->C / options()->rho_0;
    // calcul de la temperature en fonction de la chaleur specifique
    m_temperature[ev] = internal_energy / specific_heat;
    m_dpde[ev] = 0.;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real MieGruneisenEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real MieGruneisenEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real MieGruneisenEOSService::getSpecificHeatCst(IMeshEnvironment* env) { return options()->specificHeat();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_MIEGRUNEISENEOS(MieGruneisen, MieGruneisenEOSService);
