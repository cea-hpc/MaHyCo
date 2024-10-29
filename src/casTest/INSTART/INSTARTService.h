// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef INSTARTSERVICE_H
#define INSTARTSERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/INSTART/INSTART_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * class liée au cas test de INSTART
 */
class INSTARTService 
: public ArcaneINSTARTObject
{
public:
  /** Constructeur de la classe */
  INSTARTService(const ServiceBuildInfo & sbi)
    : ArcaneINSTARTObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~INSTARTService() {};

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
