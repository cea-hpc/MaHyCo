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
#define CINETEST22_COUPLAGE_SINIT  __sinit_MOD_cinetest22_couplage_sinit
#define S_CALC_CINE_VE __scalc_MOD_s_calc_cine_ve
    
    // lecture du fichier de coef
extern void  S_LEC_COEF(char* fich);
extern void  CINETEST22_COUPLAGE_SINIT(double *rho_std, double *ene_std);
extern void  S_CALC_CINE_VE(Integer *nbmail,
                double *dtime,
                double *rho,
                double *ene,
                double *Pres,
                double *Temp,
                double *Frac1,
                double *Frac2,
                double *Frac3,
                double *Frac4,
                double *Frac5,
                double *Frac6,
                double *dpde,
                double *cv,
                double *cs2,
                double *conv);
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
   *  Initialise de certains elements de l'équation d'état en reprise. 
   */
  virtual void ReinitEOS(IMeshEnvironment* env);
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
   *  Applique un endommagement dans la maille (pression nulle)
   *  si la maille est endessous de la tension limite. 
   */
  virtual void Endommagement(IMeshEnvironment* env);
  /** 
   *  Renvoie la constante adiabatic de l'environnement. 
   */
  virtual Real getAdiabaticCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la constante tension limit de l'environnement. 
   */
  virtual Real getTensionLimitCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la chaleur spécifique de l'environnement. 
   */
  virtual Real getSpecificHeatCst(IMeshEnvironment* env);
  /** 
   *  Renvoie le rapport seuil de densité  
   */
  virtual Real getdensityDamageThresold(IMeshEnvironment* env);
private:
  // taille des tableaux = nombre de maille de l'environement
  Integer nbmail;
  double *m_dtime,*m_rho,*m_ene,*m_Pres,*m_Temp;
  double *m_Frac1,*m_Frac2,*m_Frac3,*m_Frac4,*m_Frac5,*m_Frac6;
  double *m_pente_dpde,*m_Cv, *m_cs2,*m_conv;

};

#endif
