// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef IELASTO_H
#define IELASTO_H

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
class IElasto
{
public:
  /** Constructeur de la classe */
  IElasto() {};
  /** Destructeur de la classe */
  virtual ~IElasto() {};
  
public:
  /** 
   *  Initialise de l'elasto-plasticité au groupe de mailles passé en argument
   */
  virtual void initElasto(IMeshEnvironment* env) = 0;
  /** 
   *  Calcul des gradients de vitesse
   */
  virtual void ComputeVelocityGradient(Real delta_t) = 0;
    /** 
   *  Calcul du tenseur de déformation et de rotation
   */
  virtual void ComputeDeformationAndRotation() = 0;
  /** 
   *  Calcul du deviateur elasto-plastique au groupe de mailles passé en argument
   */
  virtual void ComputeElasticity(IMeshEnvironment* env, Real delta_t, Integer dim) = 0;
  /** 
   *  Calcul de la plasticité au groupe de mailles passé en argument
   */
  virtual void ComputePlasticity(IMeshEnvironment* env, Real delta_t, Integer dim) = 0;
  /** 
   *  Calcul du travail elasto-plastique
   */
  virtual void ComputeElastoEnergie(IMeshEnvironment* env, Real delta_t) = 0;
};

#endif
