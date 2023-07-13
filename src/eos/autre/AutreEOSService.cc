// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "AutreEOSService.h"

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
  char* fichier="ee.CineTest22#.Sn.00#.coeff";
  S_LEC_COEF(fichier);
   
  // initialisation rho_std,ene_std
  double rho_std(0.),ene_std(0.);
  //pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  
  CINETEST22_COUPLAGE_SINIT(&rho_std,&ene_std);
  
  //pinfo() << "Etat Std :" ;
  //pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  //pinfo() << "Nombre de maille pour allocations " << nbmail;
  
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
    m_cs2   = (double *)malloc(sizeof(double) * nbmail);
    m_conv  = (double *)malloc(sizeof(double) * nbmail);
    //pinfo() << " Allocations terminées : " << nbmail;
  }
  Integer imail = 0;
  // copie dans des tableaux Theia
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();

    m_dtime[imail] = m_global_deltat() * CONVERSION_DT ;
    m_rho[imail] = m_density[ev] * CONVERSION_DENSITE;
    m_ene[imail] = m_internal_energy[ev] * CONVERSION_ENERGIE;
    // sortie m_Pres[i] et m_Temp[i];
    m_Temp[imail] = m_temperature[ev] * CONVERSION_TEMPERATURE;
    m_Pres[imail] = m_pressure[ev] * CONVERSION_PRESSION ;
    m_Frac1[imail] = m_frac_phase1[ev];
    m_Frac2[imail] = m_frac_phase2[ev];
    m_Frac3[imail] = m_frac_phase3[ev];
    m_Frac4[imail] = m_frac_phase4[ev];
    m_Frac5[imail] = m_frac_phase5[ev];
    m_Frac6[imail] = m_frac_phase6[ev];
    // sortie m_dpde,m_cs2,m_conv;
    imail++;
  }
  //pinfo() << " Appel à la fonction " ;
  //pinfo() << " Valeurs envoyé à l'EOS ";
  //pinfo() << "  rho = " << m_rho[0];
  //pinfo() << "  ene = " << m_ene[0];
  //pinfo() << "  pres = " << m_Pres[0];
  //pinfo() << "  Temp = " << m_Temp[0];
  
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
    m_sound_speed[ev] = sqrt(m_cs2[imail] / CONVERSION_VITSON);
    m_temperature[ev] = m_Temp[imail] / CONVERSION_TEMPERATURE;
    m_pressure[ev] = m_Pres[imail] / CONVERSION_PRESSION;
    m_internal_energy[ev] = m_ene[imail] / CONVERSION_ENERGIE;
    imail++;
  }
  //pinfo() << " FIN DE  la fonction " ;
  //pinfo() << " Valeurs renvoyé par l'EOS ";
  //pinfo() << "  rho = " << m_rho[0];
  //pinfo() << "  ene = " << m_ene[0];
  //pinfo() << "  pres = " << m_Pres[0];
  //pinfo() << "  Temp = " << m_Temp[0];
  //pinfo() << "  Vitson = " << m_cs2[0];
  /* ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
   pinfo() << "  rho = " << m_density[ev] ;
   pinfo() << "  ene = " << m_internal_energy[ev] ;
   pinfo() << "  pres = " << m_pressure[ev] ;
   pinfo() << "  Temp = " << m_temperature[ev];
   pinfo() << "  Vitson = " << m_sound_speed[ev];
  }*/
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::ReinitEOS(IMeshEnvironment* env)
{
  char* fichier="ee.CineTest22#.Sn.00#.coeff";
  
  S_LEC_COEF(fichier);
   
  // initialisation rho_std,ene_std
  double rho_std(0.),ene_std(0.);
  //  pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  
  
  CINETEST22_COUPLAGE_SINIT(&rho_std,&ene_std);
  
  
  //  pinfo() << "Etat Std :" ;
  //  pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  
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
    m_cs2   = (double *)malloc(sizeof(double) * nbmail);
    m_conv  = (double *)malloc(sizeof(double) * nbmail);
    
  Integer ip(0),imail(0);
  // copie dans des tableaux Theia
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();

    m_dtime[imail] = m_global_deltat() * CONVERSION_DT ;
    m_rho[imail] = m_density[ev] * CONVERSION_DENSITE;
    m_ene[imail] = m_internal_energy[ev] * CONVERSION_ENERGIE;
    // sortie m_Pres[i] et m_Temp[i];
    m_Temp[imail] = m_temperature[ev] * CONVERSION_TEMPERATURE;
    m_Pres[imail] = m_pressure[ev] * CONVERSION_PRESSION ;
    m_Frac1[imail] = m_frac_phase1[ev];
    m_Frac2[imail] = m_frac_phase2[ev];
    m_Frac3[imail] = m_frac_phase3[ev];
    m_Frac4[imail] = m_frac_phase4[ev];
    m_Frac5[imail] = m_frac_phase5[ev];
    m_Frac6[imail] = m_frac_phase6[ev];
    // sortie m_dpde,m_cs2,m_conv;
    imail++;
  }
  pinfo() << " Appel à la fonction " << " nombre de maille " << imail ;
  pinfo() << " Valeurs envoyé à l'EOS ";
  pinfo() << "  rho = " << m_rho[ip];
  pinfo() << "  ene = " << m_ene[ip];
  pinfo() << "  pres = " << m_Pres[ip];
  pinfo() << "  Temp = " << m_Temp[ip];
  pinfo() << "   DT " << m_dtime[ip]; 
  pinfo() << "  Frac1 = " << m_Frac1[ip];
  pinfo() << "  Frac2 = " << m_Frac2[ip];
  pinfo() << "  Frac3 = " << m_Frac3[ip];
  pinfo() << "  Frac4 = " << m_Frac4[ip];
  pinfo() << "  Frac5 = " << m_Frac5[ip];
  
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
    m_sound_speed[ev] = sqrt(m_cs2[imail] / CONVERSION_VITSON);
    m_temperature[ev] = m_Temp[imail] / CONVERSION_TEMPERATURE;
    m_pressure[ev] = m_Pres[imail] / CONVERSION_PRESSION;
    m_internal_energy[ev] = m_ene[imail] / CONVERSION_ENERGIE;
    imail++;
  }
  //pinfo() << " FIN DE  la fonction " ;
  //pinfo() << " Valeurs renvoyé par^ l'EOS ";
  //pinfo() << "  rho = " << m_rho[ip];
  //pinfo() << "  ene = " << m_ene[ip];
  //pinfo() << "  pres = " << m_Pres[ip];
  //pinfo() << "  Temp = " << m_Temp[ip];
  //pinfo() << "  Vitson = " << m_cs2[ip];
  /* ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
   pinfo() << "  rho = " << m_density[ev] ;
   pinfo() << "  ene = " << m_internal_energy[ev] ;
   pinfo() << "  pres = " << m_pressure[ev] ;
   pinfo() << "  Temp = " << m_temperature[ev];
   pinfo() << "  Vitson = " << m_sound_speed[ev];
  }
  */
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::applyEOS(IMeshEnvironment* env)
{
  Integer ip(0), imail(0);
  // copie dans des tableaux Theia
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();

    m_dtime[imail] = m_global_deltat() * CONVERSION_DT ;
    m_rho[imail] = m_density[ev] * CONVERSION_DENSITE;
    m_ene[imail] = m_internal_energy[ev] * CONVERSION_ENERGIE;
    // sortie m_Pres[i] et m_Temp[i];
    m_Temp[imail] = m_temperature[ev] * CONVERSION_TEMPERATURE;
    m_Pres[imail] = m_pressure[ev] * CONVERSION_PRESSION ;
    m_Frac1[imail] = m_frac_phase1[ev];
    m_Frac2[imail] = m_frac_phase2[ev];
    m_Frac3[imail] = m_frac_phase3[ev];
    m_Frac4[imail] = m_frac_phase4[ev];
    m_Frac5[imail] = m_frac_phase5[ev];
    m_Frac6[imail] = m_frac_phase6[ev];
    m_conv[imail] = 0.;
    if (m_Temp[imail] < 100.) ip = imail;
    // sortie m_dpde,m_cs2,m_conv;
    imail++;
  }
  // pinfo() << " Appel à la fonction " ;
  // pinfo() << " Valeurs envoyé à l'EOS ";
  // pinfo() << "  rho = " << m_rho[ip];
  // pinfo() << "  ene = " << m_ene[ip];
  // pinfo() << "  pres = " << m_Pres[ip];
  // pinfo() << "  Temp = " << m_Temp[ip];
  
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
    m_sound_speed[ev] = sqrt(m_cs2[imail] / CONVERSION_VITSON);
    m_temperature[ev] = m_Temp[imail] / CONVERSION_TEMPERATURE;
    m_pressure[ev] = m_Pres[imail] / CONVERSION_PRESSION;
    imail++;
  }
  // pinfo() << " FIN DE  la fonction " ;
  // pinfo() << " Valeurs renvoyé par l'EOS pour la maille " << ip;
  // pinfo() << "  rho = " << m_rho[ip];
  // pinfo() << "  ene = " << m_ene[ip];
  // pinfo() << "  pres = " << m_Pres[ip];
  // pinfo() << "  Temp = " << m_Temp[ip];
  // pinfo() << "  Vitson = " << m_cs2[ip];
  // pinfo() << " Conversion vers les envcell ";
  
 
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  nbmail =1;
  m_dtime[0] = m_global_deltat()* CONVERSION_DT ;
  m_rho[0] = m_density[ev] * CONVERSION_DENSITE;
  m_ene[0] = m_internal_energy[ev] * CONVERSION_ENERGIE; 
  // sortie m_Pres[i] et m_Temp[i];
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
  m_sound_speed[ev] = sqrt(m_cs2[0] / CONVERSION_VITSON) ;
  m_temperature[ev] = m_Temp[0] / CONVERSION_TEMPERATURE;
  m_pressure[ev] = m_Pres[0] / CONVERSION_PRESSION;
  
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real AutreEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real AutreEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real AutreEOSService::getSpecificHeatCst(IMeshEnvironment* env) { return options()->specificHeat();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_AUTREEOS(Autre, AutreEOSService);
