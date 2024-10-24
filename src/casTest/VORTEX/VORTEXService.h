// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef VORTEXSERVICE_H
#define VORTEXSERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/VORTEX/VORTEX_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * class liée au cas test de Rider
 */
class VORTEXService 
: public ArcaneVORTEXObject
{
public:
  /** Constructeur de la classe */
  VORTEXService(const ServiceBuildInfo & sbi)
    : ArcaneVORTEXObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~VORTEXService() {};

public:
  virtual void initMatMono(Integer dim);
  virtual void initVarMono(  Integer dim, 
      SharedArray<double> densite_initiale, 
      SharedArray<double> energie_initiale, 
      SharedArray<double> pression_initiale, 
      SharedArray<double> temperature_initiale, 
      SharedArray<Real3> vitesse_initiale);
  virtual void initMat(Integer dim);
  virtual void initVar(  Integer dim, 
      SharedArray<double> densite_initiale, 
      SharedArray<double> energie_initiale, 
      SharedArray<double> pression_initiale, 
      SharedArray<double> temperature_initiale, 
      SharedArray<Real3> vitesse_initiale);
  virtual bool hasReverseOption();
  virtual Real getReverseParameter();
  virtual bool isInternalModel();
  virtual void initUtilisateur(Real3 vitesse_initiale);
};

#endif
