// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef OTHERSERVICE_H
#define OTHERSERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/OTHER/OTHER_axl.h"
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
#include "cartesian/interface/ICartesianMesh.h"
using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Clase liée au cas test unitaire
 */
class OTHERService 
: public ArcaneOTHERObject
{
public:
  /** Constructeur de la classe */
  OTHERService(const ServiceBuildInfo & sbi)
    : ArcaneOTHERObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~OTHERService() {};

public:
  virtual void initMatMono(Integer dim);
  virtual void initVarMono(Integer dim);
  virtual void initMat(Integer dim);
  virtual void initVar(Integer dim);
  virtual bool hasReverseOption();
  virtual Real getReverseParameter();
};

#endif
