// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef IEQUATIONOFSTATE_H
#define IEQUATIONOFSTATE_H

#include <arcane/ItemTypes.h>

#include "arcane/materials/IMeshMaterialMng.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/IMeshBlock.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/MeshEnvironmentVariableRef.h"
#include "arcane/materials/MaterialVariableBuildInfo.h"
#include "arcane/materials/CellToAllEnvCellConverter.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Interface du service du modèle de calcul de l'équation d'état.
 */
class IEquationOfState
{
public:
  /** Constructeur de la classe */
  IEquationOfState() {};
  /** Destructeur de la classe */
  virtual ~IEquationOfState() {};
  
public:
  /** 
   *  Initialise l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et l'énergie interne. 
   */
  virtual void initEOS(IMeshEnvironment* env) = 0;
   /** 
   *  Initialise de certains elements de l'équation d'état en reprise. 
   */
  virtual void ReinitEOS(IMeshEnvironment* env)  = 0;
  /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression. 
   */
  virtual void applyEOS(IMeshEnvironment* env) = 0;
    /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression pour une cellule
   */
  virtual void applyOneCellEOS(IMeshEnvironment* env, EnvCell ev) = 0;
  /** 
   *  Applique un endommagement dans la maille (pression nulle)
   *  si la maille est endessous de la tension limite. 
   */
  virtual void Endommagement(IMeshEnvironment* env) = 0;
  /** 
   *  Renvoie la constante adiabatic de l'environnement. 
   */
  virtual Real getAdiabaticCst(IMeshEnvironment* env) = 0;
  /** 
   *  Renvoie la constante tension limit de l'environnement. 
   */
  virtual Real getTensionLimitCst(IMeshEnvironment* env) = 0;
  /** 
   *  Renvoie la chaleur spécifique de l'environnement. 
   */
  virtual Real getSpecificHeatCst(IMeshEnvironment* env) = 0;
  /** 
   *  Renvoie le rapport seuil de densité  
   */
  virtual Real getdensityDamageThresold(IMeshEnvironment* env) = 0;
};

#endif
