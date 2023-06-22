// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef ONECELLSERVICE_H
#define ONECELLSERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/ONECELL/ONECELL_axl.h"
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
#include "arcane/cea/ICartesianMesh.h"
using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Clase liée au cas test unitaire
 */
class ONECELLService 
: public ArcaneONECELLObject
{
public:
  /** Constructeur de la classe */
  ONECELLService(const ServiceBuildInfo & sbi)
    : ArcaneONECELLObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~ONECELLService() {};

public:
  virtual void initMatMono(Integer dim);
  virtual void initVarMono(Integer dim, Real3 densite_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale);
  virtual void initMat(Integer dim);
  virtual void initVar(Integer dim, Real3 densite_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale);
  virtual bool hasReverseOption();
  virtual Real getReverseParameter();
};

#endif
