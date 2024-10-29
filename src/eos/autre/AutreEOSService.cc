// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#include "AutreEOSService.h"

#include <string> 

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#define CONVERSION_DT 1.e3           // (de s à ms) 
#define CONVERSION_DENSITE 1.e-3   // (de kg/m3 à g/cm3)
#define CONVERSION_PRESSION 1.e-9  // (de Pa à GPa)
#define CONVERSION_ENERGIE  1.e-6  // (de J/Kg à MJ/kg)
#define CONVERSION_TEMPERATURE 1.  // (de K à K)
#define CONVERSION_VITSON 1.e-6       // (de m2/s2 à m2.ms-2)

void AutreEOSService::initEOS(IMeshEnvironment* env)
{
  bool impression = false;
  String fichier_string = options()->fichierCoeff();
  const char* fichier= fichier_string.localstr();
  pinfo() << "Debut de S_LEC_COEF";
  
  S_LEC_COEF((char*) fichier);
   
  // initialisation rho_std,ene_std
  double rho_std(0.),ene_std(0.);
  if (impression) 
    pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  
  CINETEST22_COUPLAGE_SINIT(&rho_std,&ene_std);
  
  
  if (nbmail == 0) {
    // premier et seul passage 
    // au lieu taille des tableaux = nombre de maille de l'environement
    // on fait au plus simple
    // taille des tableaux = nombre de maille du calcul
    nbmail = allCells().size();
    //pinfo() << "Nombre de maille pour allocations " << nbmail;
    m_dtime = (double *)malloc(sizeof(double) * nbmail);
    m_rho   = (double *)malloc(sizeof(double) * nbmail);
    m_ene   = (double *)malloc(sizeof(double) * nbmail);
    m_Pres  = (double *)malloc(sizeof(double) * nbmail);
    m_Temp  = (double *)malloc(sizeof(double) * nbmail);
    m_Frac1 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac2 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac3 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac4 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac5 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac6 = (double *)malloc(sizeof(double) * nbmail);
    m_pente_dpde  = (double *)malloc(sizeof(double) * nbmail);
    m_Cv    = (double *)malloc(sizeof(double) * nbmail);
    m_cs2   = (double *)malloc(sizeof(double) * nbmail);
    m_conv  = (double *)malloc(sizeof(double) * nbmail);
  }
  Integer ip(0),imail(0);
  // copie dans des tableaux Theia
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();

    m_dtime[imail] = m_global_deltat() * CONVERSION_DT ;
    // solution esther : 
    m_dtime[imail] = 1. * CONVERSION_DT ;
    m_rho[imail] = m_density[ev] * CONVERSION_DENSITE;
    m_ene[imail] = m_internal_energy[ev] * CONVERSION_ENERGIE;
    // sortie m_Pres[i] et m_Temp[i];
    m_Temp[imail] = m_temperature[ev] * CONVERSION_TEMPERATURE;
    m_Pres[imail] = m_pressure[ev] * CONVERSION_PRESSION;
    m_Frac1[imail] = m_frac_phase1[ev];
    m_Frac2[imail] = m_frac_phase2[ev];
    m_Frac3[imail] = m_frac_phase3[ev];
    m_Frac4[imail] = m_frac_phase4[ev];
    m_Frac5[imail] = m_frac_phase5[ev];
    m_Frac6[imail] = m_frac_phase6[ev];
    m_Cv[imail] = m_cv[ev] * CONVERSION_ENERGIE/CONVERSION_TEMPERATURE;
    // sortie m_dpde,m_cs2,m_conv;
    imail++;
  }
  if (impression) { 
    pinfo() << " Appel à la fonction " ;
    pinfo() << " Valeurs envoyé à l'EOS ";
    pinfo() << "  rho = " << m_rho[0];
    pinfo() << "  ene = " << m_ene[0];
    pinfo() << "  pres = " << m_Pres[0];
    pinfo() << "  Temp = " << m_Temp[0];
    pinfo() << "  Frac_phase = " << m_Frac1[0] << " , " << m_Frac2[0] << " , " << m_Frac3[0];
  }
  
  S_CALC_CINE_VE( &nbmail,
                 m_dtime,
                 m_rho,
                 m_ene,
                 m_Pres,
                 m_Temp,
                 m_Frac1,
                 m_Frac2,
                 m_Frac3,
                 m_Frac4,
                 m_Frac5,
                 m_Frac6,
                 m_pente_dpde,
                 m_Cv,
                 m_cs2,
                 m_conv);
  imail = 0;
  // Retour Pression et vitesse du son
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    m_frac_phase1[ev] = m_Frac1[imail];
    m_frac_phase2[ev] = m_Frac2[imail];
    m_frac_phase3[ev] = m_Frac3[imail];
    m_frac_phase4[ev] = m_Frac4[imail];
    m_frac_phase5[ev] = m_Frac5[imail];
    m_frac_phase6[ev] = m_Frac6[imail];
    m_dpde[ev] = m_pente_dpde[imail] / (CONVERSION_PRESSION/CONVERSION_ENERGIE);
    m_cv[ev] = m_Cv[imail] / (CONVERSION_ENERGIE/CONVERSION_TEMPERATURE);
    m_sound_speed[ev] = sqrt(m_cs2[imail] / CONVERSION_VITSON);
    m_temperature[ev] = m_Temp[imail] / CONVERSION_TEMPERATURE;
    m_pressure[ev] = m_Pres[imail] / CONVERSION_PRESSION;
    m_internal_energy[ev] = m_ene[imail] / CONVERSION_ENERGIE;
    m_density_0[ev] = m_density[ev];
    m_internal_energy_0[ev] =  m_internal_energy[ev];
    imail++;
  }
  
  if (impression) {
    pinfo() << " FIN de l'init " ;
    pinfo() << " -------------- " ;
    pinfo() << " Valeurs renvoyé par l'EOS ";
    pinfo() << "  Densité          = " << m_rho[0];
    pinfo() << "  Energie Interne  = " << m_ene[0];
    pinfo() << "  Pression         = " << m_Pres[0];
    pinfo() << "  Température      = " << m_Temp[0];
    pinfo() << "  Vitesse du son   = " << m_cs2[0];
    pinfo() << "  Fraction de la phase 1 (Beta)    = " << m_Frac1[0];
    pinfo() << "  Fraction de la phase 2 (Gamma)   = " << m_Frac2[0];
    pinfo() << "  Fraction de la phase 3 (Liquide) = " << m_Frac3[0];
    pinfo() << " ------------------------------------------------- " ;
  }
    
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::ReinitEOS(IMeshEnvironment* env)
{
  bool impression = false;
  String fichier_string = options()->fichierCoeff();
  const char* fichier= fichier_string.localstr();;
  S_LEC_COEF((char*) fichier);
   
  // initialisation rho_std,ene_std
  double rho_std(0.),ene_std(0.);
  if (impression) 
    pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  
  
  CINETEST22_COUPLAGE_SINIT(&rho_std,&ene_std);
  
  // au lieu taille des tableaux = nombre de maille de l'environement
  // on fait au plus simple
  // taille des tableaux = nombre de maille du calcul
    nbmail = allCells().size();
    m_dtime = (double *)malloc(sizeof(double) * nbmail);
    m_rho   = (double *)malloc(sizeof(double) * nbmail);
    m_ene   = (double *)malloc(sizeof(double) * nbmail);
    m_Pres  = (double *)malloc(sizeof(double) * nbmail);
    m_Temp  = (double *)malloc(sizeof(double) * nbmail);
    m_Frac1 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac2 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac3 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac4 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac5 = (double *)malloc(sizeof(double) * nbmail);
    m_Frac6 = (double *)malloc(sizeof(double) * nbmail);
    m_pente_dpde  = (double *)malloc(sizeof(double) * nbmail);
    m_Cv    = (double *)malloc(sizeof(double) * nbmail);
    m_cs2   = (double *)malloc(sizeof(double) * nbmail);
    m_conv  = (double *)malloc(sizeof(double) * nbmail);
    
  Integer ip(0),imail(0);
  // copie dans des tableaux Theia
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    m_dtime[imail] = m_global_deltat() * CONVERSION_DT ;
    m_rho[imail] = m_density[ev] * CONVERSION_DENSITE;
    m_ene[imail] = m_internal_energy[ev] * CONVERSION_ENERGIE;
    // sortie m_Pres[i] et m_Temp[i];
    m_Temp[imail] = m_temperature[ev] * CONVERSION_TEMPERATURE;
    m_Pres[imail] = m_pressure[ev] * CONVERSION_PRESSION;
    m_Frac1[imail] = m_frac_phase1[ev];
    m_Frac2[imail] = m_frac_phase2[ev];
    m_Frac3[imail] = m_frac_phase3[ev];
    m_Frac4[imail] = m_frac_phase4[ev];
    m_Frac5[imail] = m_frac_phase5[ev];
    m_Frac6[imail] = m_frac_phase6[ev];
    m_Cv[imail] = m_cv[ev] * CONVERSION_ENERGIE/CONVERSION_TEMPERATURE;
    // sortie m_dpde,m_cs2,m_conv;
    imail++;
  }
  
  if (impression) { 
    pinfo() << " FIN de la Re-init " ;
    pinfo() << " -------------- " ;
    pinfo() << " Valeurs renvoyé par l'EOS ";
    pinfo() << "  Densité          = " << m_rho[0];
    pinfo() << "  Energie Interne  = " << m_ene[0];
    pinfo() << "  Pression         = " << m_Pres[0];
    pinfo() << "  Température      = " << m_Temp[0];
    pinfo() << "  Vitesse du son   = " << m_cs2[0];
    pinfo() << "  Fraction de la phase 1 (Beta)    = " << m_Frac1[0];
    pinfo() << "  Fraction de la phase 2 (Gamma)   = " << m_Frac2[0];
    pinfo() << "  Fraction de la phase 3 (Liquide) = " << m_Frac3[0];
    pinfo() << " ------------------------------------------------- " ;
  }
  
  S_CALC_CINE_VE( &nbmail,
                 m_dtime,
                 m_rho,
                 m_ene,
                 m_Pres,
                 m_Temp,
                 m_Frac1,
                 m_Frac2,
                 m_Frac3,
                 m_Frac4,
                 m_Frac5,
                 m_Frac6,
                 m_pente_dpde,
                 m_Cv,
                 m_cs2,
                 m_conv);
  
  imail = 0;
  // Retour Pression et vitesse du son
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    m_frac_phase1[ev] = m_Frac1[imail];
    m_frac_phase2[ev] = m_Frac2[imail];
    m_frac_phase3[ev] = m_Frac3[imail];
    m_frac_phase4[ev] = m_Frac4[imail];
    m_frac_phase5[ev] = m_Frac5[imail];
    m_frac_phase6[ev] = m_Frac6[imail];
    m_dpde[ev] = m_pente_dpde[imail] / (CONVERSION_PRESSION/CONVERSION_ENERGIE);
    m_cv[ev] = m_Cv[imail] / (CONVERSION_ENERGIE/CONVERSION_TEMPERATURE);
    m_sound_speed[ev] = sqrt(m_cs2[imail] / CONVERSION_VITSON);
    m_temperature[ev] = m_Temp[imail] / CONVERSION_TEMPERATURE;
    m_pressure[ev] = m_Pres[imail] / CONVERSION_PRESSION;
    m_internal_energy[ev] = m_ene[imail] / CONVERSION_ENERGIE;
    imail++;
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::applyEOS(IMeshEnvironment* env)
{
  bool impression = false;
  Integer ip(0), imail(0);
  // copie dans des tableaux Theia
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();
    if (m_maille_endo[ev.globalCell()] == 0) {
        m_dtime[imail] = m_global_deltat() * CONVERSION_DT ;
        m_rho[imail] = m_density[ev] * CONVERSION_DENSITE;
        m_ene[imail] = m_internal_energy[ev] * CONVERSION_ENERGIE;
        // sortie m_Pres[imail] et m_Temp[imail]  mais on les envoie à S_CALC_CINE_VE
        // pour initialiser le newton plus proprement
        m_Temp[imail] = m_temperature[ev] * CONVERSION_TEMPERATURE;
        // m_Pres[imail] = std::max(10. , m_pressure[ev] * CONVERSION_PRESSION) ;pour les forts Bij voir avec Greg
        m_Pres[imail] = m_pressure[ev] * CONVERSION_PRESSION ;
        m_Frac1[imail] = m_frac_phase1[ev]; 
        m_Frac2[imail] = m_frac_phase2[ev];
        m_Frac3[imail] = m_frac_phase3[ev];
        m_Frac4[imail] = m_frac_phase4[ev];
        m_Frac5[imail] = m_frac_phase5[ev];
        m_Frac6[imail] = m_frac_phase6[ev];
        m_Cv[imail] = m_cv[ev] * CONVERSION_ENERGIE/CONVERSION_TEMPERATURE;
        m_conv[imail] = 0.;
        if (m_Temp[imail] < 100.) ip = imail;
        // sortie m_dpde,m_cs2,m_conv;
        imail++;
    }
  }
  if (impression) { 
    pinfo() << " Appel à la fonction " ;
    pinfo() << " Valeurs envoyé à l'EOS ";
    pinfo() << " Frac 1 = " << m_Frac1[ip];
    pinfo() << " Frac 2 = " << m_Frac2[ip];
    pinfo() << " Frac 3 = " << m_Frac3[ip];
    pinfo() << "  Energie Interne  = " << m_ene[ip];
    pinfo() << "  Pression         = " << m_Pres[ip];
    pinfo() << "  Température      = " << m_Temp[ip]; 
    pinfo() << "  Densité          = "  << m_rho[ip];
    pinfo() << "  Pas de Temps     = " << m_dtime[ip];
    pinfo() << "  ------------------------------- " ;
  }
  
  S_CALC_CINE_VE( &nbmail,
                 m_dtime,
                 m_rho,
                 m_ene,
                 m_Pres,
                 m_Temp,
                 m_Frac1,
                 m_Frac2,
                 m_Frac3,
                 m_Frac4,
                 m_Frac5,
                 m_Frac6,
                 m_pente_dpde,
                 m_Cv,
                 m_cs2,
                 m_conv);
  
  
  if (impression) { 
    pinfo() << " Valeurs renvoyé à l'EOS ";
    pinfo() << " Frac 1 = " << m_Frac1[ip];
    pinfo() << " Frac 2 = " << m_Frac2[ip];
    pinfo() << " Frac 3 = " << m_Frac3[ip];
    pinfo() << "  Energie Interne  = " << m_ene[ip];
    pinfo() << "  Pression         = " << m_Pres[ip];
    pinfo() << "  Température      = " << m_Temp[ip];
    pinfo() << "  Densité          = "  << m_rho[ip];
    pinfo() << "  ------------------------------- " ;
  }
  
  imail = 0;
  // Retour Pression et vitesse du son
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell; 
    if (m_maille_endo[ev.globalCell()] == 0) {
        m_frac_phase1[ev] = m_Frac1[imail];
        m_frac_phase2[ev] = m_Frac2[imail];
        m_frac_phase3[ev] = m_Frac3[imail];
        m_frac_phase4[ev] = m_Frac4[imail];
        m_frac_phase5[ev] = m_Frac5[imail];
        m_frac_phase6[ev] = m_Frac6[imail];
        m_dpde[ev] = m_pente_dpde[imail] / (CONVERSION_PRESSION/CONVERSION_ENERGIE);
        m_cv[ev] = m_Cv[imail] / (CONVERSION_ENERGIE/CONVERSION_TEMPERATURE);
        m_sound_speed[ev] = sqrt(m_cs2[imail] / CONVERSION_VITSON);
        m_temperature[ev] = m_Temp[imail] / CONVERSION_TEMPERATURE;
        m_pressure[ev] = m_Pres[imail] / CONVERSION_PRESSION;
        imail++;
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/  

void AutreEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  if (m_maille_endo[ev.globalCell()] == 1) return;
  nbmail =1;
  m_dtime[0] = m_global_deltat()* CONVERSION_DT ;
  m_rho[0] = m_density[ev] * CONVERSION_DENSITE;
  m_ene[0] = m_internal_energy[ev] * CONVERSION_ENERGIE; 
  // sortie m_Pres[0] et m_Temp[0] mais on les envoie à S_CALC_CINE_VE
  // pour initialiser le newton plus proprement
  m_Temp[0] = m_temperature[ev] * CONVERSION_TEMPERATURE;
  // m_Pres[0] = std::max(10. , m_pressure[ev] * CONVERSION_PRESSION) ; pour les forts Bij voir avec Greg
  m_Pres[0] = m_pressure[ev] * CONVERSION_PRESSION ;
  m_Frac1[0] = m_frac_phase1[ev];
  m_Frac2[0] = m_frac_phase2[ev];
  m_Frac3[0] = m_frac_phase3[ev];
  m_Frac4[0] = m_frac_phase4[ev];
  m_Frac5[0] = m_frac_phase5[ev];
  m_Frac6[0] = m_frac_phase6[ev];
  
  S_CALC_CINE_VE( &nbmail,
                 m_dtime,
                 m_rho,
                 m_ene,
                 m_Pres,
                 m_Temp,
                 m_Frac1,
                 m_Frac2,
                 m_Frac3,
                 m_Frac4,
                 m_Frac5,
                 m_Frac6,
                 m_pente_dpde,
                 m_Cv,
                 m_cs2,
                 m_conv);
  
  // Retour Pression et vitesse du son
  m_frac_phase1[ev] = m_Frac1[0];
  m_frac_phase2[ev] = m_Frac2[0];
  m_frac_phase3[ev] = m_Frac3[0];
  m_frac_phase4[ev] = m_Frac4[0];
  m_frac_phase5[ev] = m_Frac5[0];
  m_frac_phase6[ev] = m_Frac6[0];
  m_dpde[ev] = m_pente_dpde[0];
  m_cv[ev] = m_Cv[0] / (CONVERSION_ENERGIE/CONVERSION_TEMPERATURE);
  m_sound_speed[ev] = sqrt(m_cs2[0] / CONVERSION_VITSON) ;
  m_temperature[ev] = m_Temp[0] / CONVERSION_TEMPERATURE;
  m_pressure[ev] = m_Pres[0] / CONVERSION_PRESSION;
  
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::Endommagement(IMeshEnvironment* env)
{
  Real damage_thresold = options()->tensionDamageThresold();
  Real density_thresold = options()->densityDamageThresold();
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;  
    Cell cell = ev.globalCell(); 
    if (m_maille_endo[ev.globalCell()] == 0) {
        // Maille saine : verification des seuils 
        if (m_pressure[ev] < damage_thresold || m_density[cell]/m_density_0[cell] < density_thresold) {
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
Real AutreEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real AutreEOSService::getdensityDamageThresold(IMeshEnvironment* env) { return options()->densityDamageThresold();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_AUTREEOS(Autre, AutreEOSService);
