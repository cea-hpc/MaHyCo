#include "VORTEXService.h"


void VORTEXService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0;
  }
}
void VORTEXService::initMat(Integer dim)  {
    
  if (options()->casTest == Vortex)  {
    initMatMono(dim);
    return;
  } 
}
void VORTEXService::initVarMono(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale)  {
    
  Real3 Xb = {5.0, 0.0, 0.};
  Real beta = 5.;
  Real gamma = 1.4;
  Real expo, func, r;
        
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // parametres maille

    // Air partout
    // Rayon
    r = sqrt((m_cell_coord[cell][0] - Xb[0]) *
                        (m_cell_coord[cell][0] - Xb[0]) +
                    (m_cell_coord[cell][1] - Xb[1]) *
                        (m_cell_coord[cell][1] - Xb[1]));
    
    expo = 1./(gamma-1.);
    func = (1-(beta*beta*(gamma-1.)/(8*gamma*Pi*Pi))*exp(1-r*r) );
    m_density[cell] = std::pow(func, expo);
    m_pressure[cell] = m_density[cell] * func;

  }
  ENUMERATE_NODE(inode, allNodes()){ 
  
    Real3 dd  = m_node_coord[inode] -Xb;
    r = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1]);
    func = exp((1-r*r)/2.) * beta /(2*Pi);
    m_velocity[inode].x = - dd[1] * func ;
    m_velocity[inode].y =  dd[0] * func ;
 
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
    
    if (m_node_coord[inode].y < 0.) pinfo() << "vit Y " << m_velocity[inode].y ;
  }
}
void VORTEXService::initVar(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale)  {

  if (options()->casTest == Vortex)  {
    initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, temperature_initiale, vitesse_initiale);
    return;
  } 
 
}

/*---------------------------------------------------------------------------*/

bool VORTEXService::hasReverseOption() { return options()->reverseOption;}
Real VORTEXService::getReverseParameter() { return options()->parametre;}
bool VORTEXService::isInternalModel() { return true; }
void VORTEXService::initUtilisateur() {}

ARCANE_REGISTER_SERVICE_VORTEX(VORTEX, VORTEXService);
