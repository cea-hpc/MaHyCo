// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef IELASTOMU_H
#define IELASTOMU_H

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
class IElastoMu
{
public:
  /** Constructeur de la classe */
  IElastoMu() {};
  /** Destructeur de la classe */
  virtual ~IElastoMu() {};
  
public:
  /** 
   *  Renvoie la constante Mu de l'environnement. 
   */
  virtual Real getShearModulus(IMeshEnvironment* env, EnvCell ev) = 0;

   Real f_ram(Real x) {
     if (x < 0.0) return 1.;
     else if (x > 1.0) return 0.;
     else return (1.0+2.0*pow(x,3)-3.0*pow(x,2));
    };  
};

#endif
