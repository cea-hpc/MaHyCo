// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SEDOVSERVICE_H
#define SEDOVSERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/SEDOV/SEDOV_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente 
 */
class SEDOVService 
: public ArcaneSEDOVObject
{
public:
  /** Constructeur de la classe */
  SEDOVService(const ServiceBuildInfo & sbi)
    : ArcaneSEDOVObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~SEDOVService() {};

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
