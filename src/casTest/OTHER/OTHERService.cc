#include "OTHERService.h"


void OTHERService::initMatMono(Integer dim)  {
}
void OTHERService::initMat(Integer dim)  {
    
  // rayon interne et externe
  double rb(0.5);
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;

    double r = m_cell_coord[cell][0];
    m_materiau[cell] = 0;
    
    if (r < rb) {
      // maille pure autre mat
    m_materiau[cell] = 1;
    } 
  }
}
void OTHERService::initVarMono(Integer dim)  {
}
void OTHERService::initVar(Integer dim)  {

  // rayon interne et externe
  double rb(0.5);
        
  info() << " boucle sur les mailles";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // pseudo-viscosité 
    m_pseudo_viscosity[cell] = 0.;
    // parametres maille
 
    // Air partout
    m_density[cell] = 0.1;
    m_pressure[cell] = 1.;
    m_internal_energy[cell] = 1.;
    // bulle surchargera l'aire
    // centre de la bulle
    double r = m_cell_coord[cell][0];
    if (r < rb) {
      // maille pure de bulle
      m_density[cell] = 1.;
    } 
    m_fracvol[cell] = 1.;
    m_mass_fraction[cell] = 1.;
  }
  info() << " boucle sur les noeuds";
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};    
    //m_velocity[inode].x = 1.;  
    //m_velocity[inode].y = 1.;
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
  info() << " fin de boucle sur les noeuds";
}

/*---------------------------------------------------------------------------*/

bool OTHERService::hasReverseOption() { return options()->reverseOption;}
Real OTHERService::getReverseParameter() { return options()->parametre;}

ARCANE_REGISTER_SERVICE_OTHER(OTHER, OTHERService); 
