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
void ONECELLService::initVarMono(Integer dim)  {
    
  info() << " Initialisation ONECELL : boucle sur les mailles";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    m_density[cell] = 7286.62;
    m_pressure[cell] = 1.e5;
    m_internal_energy[cell] = 7.6961184899;
    // pseudo-viscosité 
    m_pseudo_viscosity[cell] = 0.;
    m_fracvol[cell] = 1.;
    m_mass_fraction[cell] = 1.;
    m_frac_phase1[cell] = 1.;
    m_frac_phase2[cell] = 0.;
    m_frac_phase3[cell] = 0.;
    m_frac_phase4[cell] = 0.;
    m_frac_phase5[cell] = 0.;
    m_frac_phase6[cell] = 0.;
    m_temperature[cell] = 300.001;
  }
    
}
void ONECELLService::initVar(Integer dim)  {
  initVarMono(dim);
}

/*---------------------------------------------------------------------------*/

bool ONECELLService::hasReverseOption() { return options()->reverseOption;}
Real ONECELLService::getReverseParameter() { return options()->parametre;}

ARCANE_REGISTER_SERVICE_ONECELL(ONECELL, ONECELLService); 
