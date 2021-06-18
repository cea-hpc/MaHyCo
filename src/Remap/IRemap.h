// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef IREMAP_H
#define IREMAP_H

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
  virtual void appliRemap(Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * resize les variables du remap
   **/
  virtual void resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * calcul des gradients aux faces
   **/
  virtual void computeGradPhiFace(Integer idir, Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * calcul des gradients aux faces ou flux aux faces 
   **/
  virtual void computeGradPhiCell(Integer idir, Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * calcul des flux aux faces des cellules 
   **/
  virtual void computeUpwindFaceQuantitiesForProjection(Integer idir, Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * calcul des valeurs aux cellules 
   **/
  virtual void computeUremap(Integer idir, Integer nb_vars_to_project, Integer nb_env) = 0;
  /**
   * synchronisation des valeurs aux cellules 
   **/
  virtual void synchronizeUremap() = 0;

  
  virtual Integer getOrdreProjection() = 0;
  virtual bool hasProjectionPenteBorne() = 0;    
  virtual bool hasConservationEnergieTotale() = 0;    
};

#endif
