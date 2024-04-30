﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef CHOCBULLESERVICE_H
#define CHOCBULLESERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/CHOCBULLE/CHOCBULLE_axl.h"
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
 * class liée au cas test de CHOCBULLE
 */
class CHOCBULLEService 
: public ArcaneCHOCBULLEObject
{
public:
  /** Constructeur de la classe */
  CHOCBULLEService(const ServiceBuildInfo & sbi)
    : ArcaneCHOCBULLEObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~CHOCBULLEService() {};

public:
  virtual void initMatMono(Integer dim);
  virtual void initVarMono(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale);
  virtual void initMat(Integer dim);
  virtual void initVar(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale);
  virtual bool hasReverseOption();
  virtual Real getReverseParameter();
  virtual bool isInternalModel();
  virtual void initUtilisateur(Real3 vitesse_initiale);
private:
  // valeur aléatoire pour les bulles
  float *m_rand1, *m_rand2, *m_rand3, *m_rand4;  
  std::vector<double> m_posx;
  std::vector<double> m_posy;
  std::vector<double> m_posz;
  std::vector<double> m_posr;
  void lectureFichier(const std::string& nomFichier);
};

#endif
