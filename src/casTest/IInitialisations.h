﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef IINITIALISATION_H
#define IINITIALISATION_H

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

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Interface du service du modèle de calcul de l'équation d'état.
 */
class IInitialisations
{
public:
  /** Constructeur de la classe */
  IInitialisations() {};
  /** Destructeur de la classe */
  virtual ~IInitialisations() {};
  
public:
  /** 
   *  Initialise les cas test en mono ou multimatériau
   */
  virtual void initMatMono(Integer dim)  = 0;
  virtual void initVarMono(  Integer dim, 
    SharedArray<double> densite_initiale, 
    SharedArray<double> energie_initiale, 
    SharedArray<double> pression_initiale, 
    SharedArray<double> temperature_initiale, 
    SharedArray<Real3> vitesse_initiale) = 0;
  virtual void initMat(Integer dim)  = 0;
  virtual void initVar(  Integer dim, 
    SharedArray<double> densite_initiale, 
    SharedArray<double> energie_initiale, 
    SharedArray<double> pression_initiale, 
    SharedArray<double> temperature_initiale, 
    SharedArray<Real3> vitesse_initiale) = 0;
  virtual bool hasReverseOption() = 0;
  virtual Real getReverseParameter() = 0;
  virtual bool isInternalModel() = 0;
  virtual void initUtilisateur(Real3 vitesse_initiale) = 0;
};

#endif
