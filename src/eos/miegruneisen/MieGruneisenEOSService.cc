// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#include "MieGruneisenEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MieGruneisenEOSService::initEOS(IMeshEnvironment* env)
{
  // Récupère le cv 
  Real cv = options()->specificHeat;
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
    m_internal_energy[ev] =  m_pressure[ev] / options()->gamma0;
    m_sound_speed[ev] = 1;
    if (eref == 0) eref = m_internal_energy[ev];
    m_temperature[ev] = (m_internal_energy[ev] - eref) / cv + tref;
    m_density_0[ev] = m_density[ev];
     m_internal_energy_0[ev] =  m_internal_energy[ev];
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
  // Récupère le cv 
  Real cv = options()->specificHeat;
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
        Real gruneisen = options()->gamma0 * J + options()->a * (1-J);
        Real numerateur  = options()->rho0 * pow(options()->c0, 2) * mu * ( 1 + ( 1. -  options()->gamma0/2.) * mu - options()->a * pow(mu, 2)/2.);
        Real denominateur = pow(( 1. - (options()->s1 -1.) * mu - options()->s2 * pow(mu, 2)/ (1+mu) - options()->s2  * pow(mu, 3)/pow(1.+mu, 2)), 2);
        m_pressure[ev] = numerateur / denominateur + (options()->gamma0 + options()->a * mu) * m_internal_energy[ev];
        // calcul de la temperature en fonction de la Capacité thermique isochore
        m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / cv + m_temperature_n[ev];
        m_dpde[ev] = (options()->gamma0 + options()->a * mu);
        // vitesse du son 
        J = options()->rho0 / (density - 1.e-7);
        mu = 1./J -1.;
        gruneisen = options()->gamma0 * J + options()->a * (1-J);
        numerateur  = options()->rho0 * pow(options()->c0, 2) * mu * ( 1 + ( 1. -  options()->gamma0/2.) * mu - options()->a * pow(mu, 2)/2.);
        denominateur = pow(( 1. - (options()->s1 -1.) * mu - options()->s2 * pow(mu, 2)/ (1+mu) - options()->s2  * pow(mu, 3)/pow(1.+mu, 2)), 2);
        Real pres1 = numerateur / denominateur + (options()->gamma0 + options()->a * mu) * m_internal_energy[ev];
        J = options()->rho0 / (density + 1.e-7);
        mu = 1./J -1.;
        gruneisen = options()->gamma0 * J + options()->a * (1-J);
        numerateur  = options()->rho0 * pow(options()->c0, 2) * mu * ( 1 + ( 1. -  options()->gamma0/2.) * mu - options()->a * pow(mu, 2)/2.);
        denominateur = pow(( 1. - (options()->s1 -1.) * mu - options()->s2 * pow(mu, 2)/ (1+mu) - options()->s2  * pow(mu, 3)/pow(1.+mu, 2)), 2);
        Real pres2 = numerateur / denominateur + (options()->gamma0 + options()->a * mu) * m_internal_energy[ev];       
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
  // Récupère le cv 
  Real cv = options()->specificHeat;
  // Calcul de la pression,la vitesse du son  et le gradient de pression pour une maille donnée
  Real internal_energy = m_internal_energy[ev];
  Real density = m_density[ev];
  if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
  Real J = options()->rho0 / density ;
    Real mu = 1./J -1.;
    Real gruneisen = options()->gamma0 * J + options()->a * (1-J);
    Real numerateur  = options()->rho0 * pow(options()->c0, 2) * mu * ( 1 + ( 1. -  options()->gamma0/2.) * mu - options()->a * pow(mu, 2)/2.);
    Real     denominateur = pow(( 1. - (options()->s1 -1.) * mu - options()->s2 * pow(mu, 2)/ (1+mu) - options()->s2  * pow(mu, 3)/pow(1.+mu, 2)), 2);
    m_pressure[ev] = numerateur / denominateur + (options()->gamma0 + options()->a * mu) * m_internal_energy[ev];
    // calcul de la temperature en fonction de la Capacité thermique isochore
    m_temperature[ev] = (m_internal_energy[ev] - m_internal_energy_n[ev]) / cv + m_temperature_n[ev];
    m_dpde[ev] = (options()->gamma0 + options()->a * mu);
    // vitesse du son 
    J = options()->rho0 / (density - 1.e-7);
    mu = 1./J -1.;
    Real pres1  = 1.;
    J = options()->rho0 / (density + 1.e-7);
    mu = 1./J -1.;
    Real pres2  = 1;     
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
Real MieGruneisenEOSService::getdensityDamageThresold(IMeshEnvironment* env) { return options()->densityDamageThresold();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
 

ARCANE_REGISTER_SERVICE_MIEGRUNEISENEOS(MieGruneisen, MieGruneisenEOSService);
