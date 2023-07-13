#include "ONECELLService.h"


void ONECELLService::initMatMono(Integer dim)  {
}
void ONECELLService::initMat(Integer dim)  {
    
  // rayon interne et externe
  double rb(0.5);
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0;
  }
}
void ONECELLService::initVarMono(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
    
  info() << " Initialisation ONECELL : boucle sur les mailles";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    m_density[cell] = densite_initiale[0];
    m_pressure[cell] = pression_initiale[0];
    m_internal_energy[cell] = energie_initiale[0];
    m_frac_phase1[cell] = 1.;
    m_frac_phase2[cell] = 0.;
    m_frac_phase3[cell] = 0.;
    m_frac_phase4[cell] = 0.;
    m_frac_phase5[cell] = 0.;
    m_frac_phase6[cell] = 0.;
  }
}
void ONECELLService::initVar(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
  initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, vitesse_initiale);
}

/*---------------------------------------------------------------------------*/

bool ONECELLService::hasReverseOption() { return options()->reverseOption;}
Real ONECELLService::getReverseParameter() { return options()->parametre;}

ARCANE_REGISTER_SERVICE_ONECELL(ONECELL, ONECELLService); 
