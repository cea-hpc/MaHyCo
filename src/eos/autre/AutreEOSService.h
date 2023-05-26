// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef AUTREEOSSERVICE_H
#define AUTREEOSSERVICE_H

#include "eos/IEquationOfState.h"
#include "eos/autre/AutreEOS_axl.h"
#include "arcane/materials/IMeshMaterialMng.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/IMeshBlock.h"
#include "arcane/materials/MeshMaterialModifier.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/MeshEnvironmentVariableRef.h"
#include "arcane/materials/MaterialVariableBuildInfo.h"
#include "arcane/materials/MeshBlockBuildInfo.h"
#include "arcane/materials/MeshEnvironmentBuildInfo.h"
#include "arcane/materials/MeshMaterialVariableDependInfo.h"
#include "arcane/materials/CellToAllEnvCellConverter.h"
#include "arcane/materials/MatCellVector.h"
#include "arcane/materials/EnvCellVector.h"
#include "arcane/materials/MatConcurrency.h"
#include "arcane/materials/MeshMaterialIndirectModifier.h"
#include "arcane/materials/MeshMaterialVariableSynchronizerList.h"
#include "arcane/materials/ComponentSimd.h"
using namespace Arcane;
using namespace Arcane::Materials;

 
#ifdef __cplusplus    
extern "C" {
#endif

#define S_LEC_COEF __m_lec_MOD_s_lec_coeff
#define TV2_CINETEST22_COUPLAGE_SINIT  __tv2_sinit_MOD_tv2_cinetest22_couplage_sinit
#define TV2_S_CALC_CINE_VE __tv2_scalc_MOD_tv2_s_calc_cine_ve
    
    // lecture du fichier de coef
extern void  S_LEC_COEF(char* fich);
extern void  TV2_CINETEST22_COUPLAGE_SINIT(float *rho_std, float *ene_std);
extern void  TV2_S_CALC_CINE_VE(Integer *nbmail,
                float *dtime,
                float *rho,
                float *ene,
                float *Pres,
                float *Temp,
                float *Frac1,
                float *Frac2,
                float *Frac3,
                float *Frac4,
                float *Frac5,
                float *Frac6,
                float *dpde,
                float *cs2,
                float *conv);
#ifdef __cplusplus       
}  
#endif

/**
 * Représente le modèle d'équation d'état <em>Gaz Parfait</em>
 */
class AutreEOSService 
: public ArcaneAutreEOSObject
{
public:
  /** Constructeur de la classe */
  AutreEOSService(const ServiceBuildInfo & sbi)
    : ArcaneAutreEOSObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~AutreEOSService() {};

public:
  /** 
   *  Initialise l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et l'énergie interne. 
   */
  virtual void initEOS(IMeshEnvironment* env);
  /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression. 
   */
  virtual void applyEOS(IMeshEnvironment* env);
   /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression pour une cellule
   */
  virtual void applyOneCellEOS(IMeshEnvironment* env, EnvCell ev);
  /** 
   *  Renvoie la constante adiabatic de l'environnement. 
   */
  virtual Real getAdiabaticCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la constante tension limit de l'environnement. 
   */
  virtual Real getTensionLimitCst(IMeshEnvironment* env);

private:
  // taille des tableaux = nombre de maille de l'environement
  Integer nbmail;
  float *m_dtime,*m_rho,*m_ene,*m_Pres,*m_Temp;
  float *m_Frac1,*m_Frac2,*m_Frac3,*m_Frac4,*m_Frac5,*m_Frac6;
  float *m_pente_dpde,*m_cs2,*m_conv;

};

#endif
