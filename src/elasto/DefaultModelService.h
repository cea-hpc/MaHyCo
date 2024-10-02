// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef DEFAULTMODELSERVICE_H
#define DEFAULTMODELSERVICE_H

#include "elasto/IElasto.h"
#include "elasto/DefaultModel_axl.h"
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
 * Représente le modèle d'élastop-plasticité
 */
class DefaultModelService 
: public ArcaneDefaultModelObject
{
public:
  /** Constructeur de la classe */
  DefaultModelService(const ServiceBuildInfo & sbi)
    : ArcaneDefaultModelObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~DefaultModelService() {};

public:
  /** 
   *  Initialise de l'elasto-plasticité au groupe de mailles passé en argument
   */
  virtual void initElasto(IMeshEnvironment* env);
  /** 
   *  Calcul des gradients de vitesse
   */
  virtual void ComputeVelocityGradient(Real delta_t);
    /** 
   *  Calcul du tenseur de déformation et de rotation
   */
  virtual void ComputeDeformationAndRotation();
  /** 
   *  Calcul du deviateur elasto-plastique au groupe de mailles passé en argument
   */
  virtual void ComputeElasticity(IMeshEnvironment* env, Real delta_t, Integer dim);
  /** 
   *  Calcul de la plasticité au groupe de mailles passé en argument
   */
  virtual void ComputePlasticity(IMeshEnvironment* env, Real delta_t, Integer dim);
  /** 
   *  Calcul du travail elasto-plastique
   */
  virtual void ComputeElastoEnergie(IMeshEnvironment* env, Real delta_t);
  /** 
   *  Renvoie la constante Mu de l'environnement. 
   */
  virtual Real getElasticCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la La limite de l'environnement. 
   */
  virtual Real getLimitElasticCst(IMeshEnvironment* env);
  
};

#endif
