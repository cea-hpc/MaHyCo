// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef EXTERNALSERVICE_H
#define EXTERNALSERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include "casTest/EXTERNAL/EXTERNAL_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Clase liée au cas test unitaire
 */
class EXTERNALService 
: public ArcaneEXTERNALObject
{
public:
  /** Constructeur de la classe */
  EXTERNALService(const ServiceBuildInfo & sbi)
    : ArcaneEXTERNALObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~EXTERNALService() {};

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
