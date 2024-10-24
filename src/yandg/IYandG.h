// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef IYANDG_H
#define IYANDG_H

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
class IYandG
{
public:
  /** Constructeur de la classe */
  IYandG() {};
  /** Destructeur de la classe */
  virtual ~IYandG() {};
  
public:
  /** 
   *  Renvoie la constante Mu de l'environnement. 
   */
  virtual Real getShearModulus(IMeshEnvironment* env, EnvCell ev) = 0;
  /** 
   *  Renvoie la La limite de l'environnement. 
   */
  virtual Real getElasticLimit(IMeshEnvironment* env, EnvCell ev) = 0;
  
};

#endif
