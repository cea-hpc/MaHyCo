// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "FictifEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FictifEOSService::initEOS(IMeshEnvironment* env)
{
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  Real eref = options()->energieRef();
  Real tref = options()->temperatureRef();
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
    if (eref == 0) eref = m_internal_energy[ev];
    m_temperature[ev] = (m_internal_energy[ev] - eref) / specific_heat + tref;
    m_density_0[ev] = m_density[ev];
    
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FictifEOSService::ReinitEOS(IMeshEnvironment* env)
{
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FictifEOSService::applyEOS(IMeshEnvironment* env)
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
    Real Tdebut(options()->tdebutPression);
    Real Tfin(options()->tfinPression);
    Real coeff(options()->valeurDependanceTemps);
    Real initvalue(options()->valeurDebutPression);
    if (m_global_time() > Tdebut ) {
         m_pressure[ev] = initvalue + coeff*m_global_time();
    }
    if (m_global_time() > Tfin ) m_pressure[ev] = 0.;
    m_sound_speed[ev] = 1. ; // pour ne pas influencer le pas de temps sqrt(adiabatic_cst * pressure / density);
    m_dpde[ev] = (adiabatic_cst - 1.) * density;
    // calcul de la temperature en fonction de la chaleur specifique
    m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / specific_heat + m_temperature_n[ev];    
    
    /* Cell cell = ev.globalCell();
    if (cell.localId() == 2000) pinfo() << env->name() << " m_temperature_n[ev] " << m_temperature_n[ev] << " m_temperature[ev] " << m_temperature[ev]; 
    if (cell.localId() == 2000) pinfo() << env->name() << " m_internal_energy_n[ev] " << m_internal_energy_n[ev] << " m_internal_energy[ev] " << m_internal_energy[ev];
    */
    // Test pas d'energie interne dan le fictif
    m_internal_energy[ev] =0.;
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FictifEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  //
  Real internal_energy = m_internal_energy[ev];
  Real density = m_density[ev];
  Real Tdebut(options()->tdebutPression);
  Real Tfin(options()->tfinPression);
  Real coeff(options()->valeurDependanceTemps);
  Real initvalue(options()->valeurDebutPression);
  if (m_global_time() > Tdebut ) {
        m_pressure[ev] = initvalue + coeff*m_global_time();
  }
  if (m_global_time() > Tfin ) m_pressure[ev] = 0.;
  m_sound_speed[ev] = 1. ; // pour ne pas influencer le pas de temps sqrt(adiabatic_cst * pressure / density);
  m_dpde[ev] = (adiabatic_cst - 1.) * density;
  // calcul de la temperature en fonction de la chaleur specifique
  m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / specific_heat + m_temperature_n[ev]; 
  // Test pas d'energie interne dan le fictif
  m_internal_energy[ev] =0.;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FictifEOSService::Endommagement(IMeshEnvironment* env)
{

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real FictifEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real FictifEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real FictifEOSService::getSpecificHeatCst(IMeshEnvironment* env) { return options()->specificHeat();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_FICTIFEOS(Fictif, FictifEOSService);
