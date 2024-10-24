// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef RIDERSERVICE_H
#define RIDERSERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/RIDER/RIDER_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * class liée au cas test de Rider
 */
class RIDERService 
: public ArcaneRIDERObject
{
public:
  /** Constructeur de la classe */
  RIDERService(const ServiceBuildInfo & sbi)
    : ArcaneRIDERObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~RIDERService() {};

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
