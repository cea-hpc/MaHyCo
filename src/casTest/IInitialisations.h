// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef IINITIALISATION_H
#define IINITIALISATION_H

#include <arcane/ItemTypes.h>

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
  virtual void initVarMono(Integer dim) = 0;
  virtual void initMat(Integer dim)  = 0;
  virtual void initVar(Integer dim) = 0;
  virtual bool hasReverseOption() = 0;
  virtual Real getReverseParameter() = 0;
};

#endif
