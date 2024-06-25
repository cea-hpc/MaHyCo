// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#include "UNIFORMService.h"

void UNIFORMService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void UNIFORMService::initVarMono(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_density[cell] = densite_initiale[0];
    m_pressure[cell] = pression_initiale[0];
    m_temperature[cell] = temperature_initiale[0];
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = vitesse_initiale[0];
  }
}

void UNIFORMService::initMat(Integer dim)  {

  // pour le moment pas de multimatériaux
  initMatMono(dim);
}

void UNIFORMService::initVar(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale)  {
    
  // pour le moment pas de multimatériaux
  initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, temperature_initiale, vitesse_initiale);
}


/*---------------------------------------------------------------------------*/


bool UNIFORMService::hasReverseOption() { return options()->reverseOption;}
Real UNIFORMService::getReverseParameter() { return options()->parametre;}
bool UNIFORMService::isInternalModel() { return true; }
void UNIFORMService::initUtilisateur(Real3 vitesse_initiale) {}

/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_UNIFORM(UNIFORM, UNIFORMService);
