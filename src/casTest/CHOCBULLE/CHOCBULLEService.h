// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef CHOCBULLESERVICE_H
#define CHOCBULLESERVICE_H

#include "TypesMahyco.h"
#include "casTest/IInitialisations.h"
#include <arcane/utils/Real3.h>

#include "casTest/CHOCBULLE/CHOCBULLE_axl.h"

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
