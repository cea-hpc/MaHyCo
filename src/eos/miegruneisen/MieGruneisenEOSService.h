// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef MIEGRUNEISENEOSSERVICE_H
#define MIEGRUNEISENEOSSERVICE_H

#include "eos/IEquationOfState.h"
#include "eos/miegruneisen/MieGruneisenEOS_axl.h"
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

/**
 * Représente le modèle d'équation d'état <em>Gaz Parfait</em>
 */
class MieGruneisenEOSService 
: public ArcaneMieGruneisenEOSObject
{
public:
  /** Constructeur de la classe */
  MieGruneisenEOSService(const ServiceBuildInfo & sbi)
    : ArcaneMieGruneisenEOSObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~MieGruneisenEOSService() {};

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
};

#endif
