﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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
  Real eref = options()->energieRef();
  Real tref = options()->temperatureRef();
  // Initialise l'énergie et la vitesse du son pour chaque maille de l'environnement
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real pressure = m_pressure[ev];
    Real density = m_density[ev];
    Real J = options()->rho0 / density ;
    Real mu = 1./J -1.;
    // Affiche l'identifiant local de la maille, la pression et la densité
    m_internal_energy[ev] = (pressure - options()->cCst * mu + options()->dCst * pow(mu,2.) + options()->sCst * pow(mu,3)) / (adiabatic_cst * density);
    m_sound_speed[ev] = math::sqrt(options()->cCst / options()->rho0);
    if (eref == 0) eref = m_internal_energy[ev];
    m_temperature[ev] = (m_internal_energy[ev] - eref) / specific_heat + tref;
    m_density_0[ev] = m_density[ev];
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MieGruneisenEOSService::ReinitEOS(IMeshEnvironment* env)
{
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
    
    if (m_maille_endo[ev.globalCell()] == 0) {
        Real internal_energy = m_internal_energy[ev];
        Real density = m_density[ev];
        if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
        Real J = options()->rho0 / density ;
        Real mu = 1./J -1.;
        m_pressure[ev] = options()->cCst * mu + options()->dCst * pow(mu,2.) + options()->sCst * pow(mu,3) +   adiabatic_cst * density * internal_energy;
        m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / specific_heat + m_temperature_n[ev];
        m_dpde[ev] = adiabatic_cst * density;  
        // vitesse du son 
        J = options()->rho0 / (density - 1.e-7);
        mu = 1./J -1.;
        Real pres1  = options()->cCst * mu + options()->dCst * pow(mu,2.) + options()->sCst * pow(mu,3) +   adiabatic_cst * density * internal_energy;
        J = options()->rho0 / (density + 1.e-7);
        mu = 1./J -1.;
        Real pres2  = options()->cCst * mu + options()->dCst * pow(mu,2.) + options()->sCst * pow(mu,3) +   adiabatic_cst * density * internal_energy;        
        if ((pres2 -  pres1) > 0.) {
            m_sound_speed[ev] = math::sqrt((pres2 -  pres1)/2.e-7);
        } else m_sound_speed[ev] = 1.e4;
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MieGruneisenEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  if (m_maille_endo[ev.globalCell()] == 1) return;
  // Récupère les constantes adiabatique et de chaleur spécifique
  Real adiabatic_cst = getAdiabaticCst(env);
  Real specific_heat = getSpecificHeatCst(env);
  // Calcul de la pression,la vitesse du son  et le gradient de pression pour une maille donnée
  Real internal_energy = m_internal_energy[ev];
  Real density = m_density[ev];
  if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
  Real J = options()->rho0 / density ;
    Real mu = 1./J -1.;
    m_pressure[ev] = options()->cCst * mu + options()->dCst * pow(mu,2.) + options()->sCst * pow(mu,3) +   adiabatic_cst * density * internal_energy;
    // calcul de la temperature en fonction de la chaleur specifique
    m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / specific_heat + m_temperature_n[ev];
    m_dpde[ev] = adiabatic_cst * density;
    // vitesse du son 
    J = options()->rho0 / (density - 1.e-7);
    mu = 1./J -1.;
    Real pres1  = options()->cCst * mu + options()->dCst * pow(mu,2.) + options()->sCst * pow(mu,3) +   adiabatic_cst * density * internal_energy;
    J = options()->rho0 / (density + 1.e-7);
    mu = 1./J -1.;
    Real pres2  = options()->cCst * mu + options()->dCst * pow(mu,2.) + options()->sCst * pow(mu,3) +   adiabatic_cst * density * internal_energy;        
    if ((pres2 -  pres1) > 0.) {
        m_sound_speed[ev] = math::sqrt((pres2 -  pres1)/2.e-7);
    } else m_sound_speed[ev] = 1.e4;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MieGruneisenEOSService::Endommagement(IMeshEnvironment* env)
{
  Real damage_thresold = options()->tensionDamageThresold();
  Real density_thresold = options()->densityDamageThresold();
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;  
    Cell cell = ev.globalCell();  
    if (m_maille_endo[ev.globalCell()] == 0) {
        // Maille saine : verification des seuils 
        if (m_pressure[cell] < damage_thresold) {
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
Real MieGruneisenEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real MieGruneisenEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real MieGruneisenEOSService::getSpecificHeatCst(IMeshEnvironment* env) { return options()->specificHeat();}
Real MieGruneisenEOSService::getdensityDamageThresold(IMeshEnvironment* env) { return options()->densityDamageThresold();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_MIEGRUNEISENEOS(MieGruneisen, MieGruneisenEOSService);
