﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef POINTTRIPLESERVICE_H
#define POINTTRIPLESERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/POINTTRIPLE/POINTTRIPLE_axl.h"
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
 * class liée au cas test du POINT TRIPLE
 */
class POINTTRIPLEService 
: public ArcanePOINTTRIPLEObject
{
public:
  /** Constructeur de la classe */
  POINTTRIPLEService(const ServiceBuildInfo & sbi)
    : ArcanePOINTTRIPLEObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~POINTTRIPLEService() {};

public:
  virtual void initMatMono(Integer dim);
  virtual void initVarMono(Integer dim, double* densite_initiale, double* energie_initiale, 
                            double* temperature_initiale, double* pression_initiale, Real3x3 vitesse_initiale);
  virtual void initMat(Integer dim);
  virtual void initVar(Integer dim, double* densite_initiale, double* energie_initiale, 
                        double* temperature_initiale, double* pression_initiale, Real3x3 vitesse_initiale);
  virtual bool hasReverseOption();
  virtual Real getReverseParameter();
  virtual bool isInternalModel();
  virtual void initUtilisateur(Real3 vitesse_initiale);
};

#endif
