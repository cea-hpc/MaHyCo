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
void OTHERService::initVarMono(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
}
void OTHERService::initVar(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {

  // rayon interne et externe
  double rb(0.5);
        
  info() << " boucle sur les mailles";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // Air partout
    m_density[cell] = densite_initiale[0];
    m_pressure[cell] = pression_initiale[0];
    m_internal_energy[cell] = energie_initiale[0];
    // bulle surchargera l'aire
    // centre de la bulle
    double r = m_cell_coord[cell][0];
    if (r < rb) {
      // maille pure de bulle
      m_density[cell] = densite_initiale[1];
    } 
  }
  info() << " boucle sur les noeuds";
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = vitesse_initiale[0];   
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
