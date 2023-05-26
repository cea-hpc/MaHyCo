// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "AutreEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
extern  void SLECCOEFSANSFICH();

#define CONVERSION_DT 1.           // (de s à s) 
#define CONVERSION_DENSITE 1.e-3   // (de kg/m3 à g/cm3)
#define CONVERSION_PRESSION 1.e-9  // (de Pa à GPa)
#define CONVERSION_ENERGIE  1.e-6  // (de J/Kg à MJ/kg)
#define CONVERSION_TEMPERATURE 1.  // (de K à K)
#define CONVERSION_VITSON 1.       // (de m/s à m/s)

void AutreEOSService::initEOS(IMeshEnvironment* env)
{
    
  // Lecture du fichier de coefficients
  S_LEC_COEF("ee.CineTest22#.Sn.00#.coeff");
   
  // initialisation rho_std,ene_std
  float rho_std(0.),ene_std(0.);
  pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  TV2_CINETEST22_COUPLAGE_SINIT(&rho_std,&ene_std);
  pinfo() << "Etat Std :" ;
  pinfo() << "rho_std=" << rho_std << "  ene_std=" << ene_std;
  pinfo() << "Nombre de maille pour allocations " << nbmail;
  
  if (nbmail == 0) {
    // premier et seul passage 
    // au lieu taille des tableaux = nombre de maille de l'environement
    // on fait au plus simple
    // taille des tableaux = nombre de maille du calcul
    nbmail = allCells().size();
    pinfo() << "Nombre de maille pour allocations " << nbmail;
    m_dtime = (float *)malloc(sizeof(float) * nbmail);
    m_rho   = (float *)malloc(sizeof(float) * nbmail);
    m_ene   = (float *)malloc(sizeof(float) * nbmail);
    m_Pres  = (float *)malloc(sizeof(float) * nbmail);
    m_Temp  = (float *)malloc(sizeof(float) * nbmail);
    m_Frac1 = (float *)malloc(sizeof(float) * nbmail);
    m_Frac2 = (float *)malloc(sizeof(float) * nbmail);
    m_Frac3 = (float *)malloc(sizeof(float) * nbmail);
    m_Frac4 = (float *)malloc(sizeof(float) * nbmail);
    m_Frac5 = (float *)malloc(sizeof(float) * nbmail);
    m_Frac6 = (float *)malloc(sizeof(float) * nbmail);
    m_pente_dpde  = (float *)malloc(sizeof(float) * nbmail);
    m_cs2   = (float *)malloc(sizeof(float) * nbmail);
    m_conv  = (float *)malloc(sizeof(float) * nbmail);
    pinfo() << " Allocations terminées : " << nbmail;
  }
  Integer imail = 0;
  // copie dans des tableaux Theia
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();

    m_dtime[imail] = m_global_deltat();
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
  pinfo() << " Apple à la fonction " ;
  TV2_S_CALC_CINE_VE( &nbmail,
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
    m_sound_speed[ev] = sqrt(m_cs2[imail]) / CONVERSION_VITSON;
    m_temperature[ev] = m_Temp[imail] / CONVERSION_TEMPERATURE;
    m_pressure[ev] = m_Pres[imail] / CONVERSION_PRESSION;
    imail++;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::applyEOS(IMeshEnvironment* env)
{
  
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AutreEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  nbmail =1;
  m_dtime[0] = m_global_deltat();
  m_rho[0] = 1./ m_density[ev];
  m_ene[0] = m_internal_energy[ev]; 
  // sortie m_Pres[i] et m_Temp[i];
  m_Frac1[0] = m_frac_phase1[ev];
  m_Frac2[0] = m_frac_phase2[ev];
  m_Frac3[0] = m_frac_phase3[ev];
  m_Frac4[0] = m_frac_phase4[ev];
  m_Frac5[0] = m_frac_phase5[ev];
  m_Frac6[0] = m_frac_phase6[ev];
  
  TV2_S_CALC_CINE_VE( &nbmail,
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
  m_sound_speed[ev] = sqrt(m_cs2[0]);
  m_temperature[ev] = m_Temp[0];
  m_pressure[ev] = m_Pres[0];
    
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real AutreEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real AutreEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_AUTREEOS(Autre, AutreEOSService);
