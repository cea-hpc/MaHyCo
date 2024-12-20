﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef IREMAP_H
#define IREMAP_H

#include <arcane/ItemTypes.h>

#include "arcane/materials/IMeshMaterialMng.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/IMeshBlock.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/MeshEnvironmentVariableRef.h"
#include "arcane/materials/MaterialVariableBuildInfo.h"
#include "arcane/materials/CellToAllEnvCellConverter.h"

#include "arcane/cea/ICartesianMesh.h"
#include "arcane/cea/CellDirectionMng.h"
#include "arcane/cea/FaceDirectionMng.h"
#include "arcane/cea/NodeDirectionMng.h"
#include "arcane/cea/CartesianConnectivity.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Interface du service du modèle de calcul de l'équation d'état.
 */
class IRemap
{
public:
  /** Constructeur de la classe */
  IRemap() {};
  /** Destructeur de la classe */
  virtual ~IRemap() {};
  
public:
    
   /**
   * main du remap
   **/
  virtual void appliRemap(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * resize les variables du remap
   **/
  virtual void resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * synchronisation des valeurs aux cellules 
   **/
  virtual void synchronizeUremap() = 0;

  
  virtual Integer getOrdreProjection() = 0;
  virtual bool hasProjectionPenteBorne() = 0; 
  virtual bool hasProjectionSimplePente() = 0;    
  virtual bool hasConservationEnergieTotale() = 0;   
  virtual String isEuler() = 0;
    /**
   * fonction final de la projection
   **/
  virtual void remapVariables(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env)  = 0;
};

#endif
