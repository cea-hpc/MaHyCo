#include "EXTERNALService.h"


void EXTERNALService::initMatMono(Integer dim)  {
}
void EXTERNALService::initMat(Integer dim)  {
   
}
void EXTERNALService::initVarMono(  Integer dim, 
  SharedArray<double> densite_initiale, 
  SharedArray<double> energie_initiale, 
  SharedArray<double> pression_initiale, 
  SharedArray<double> temperature_initiale, 
  SharedArray<Real3> vitesse_initiale)  {
}
void EXTERNALService::initVar(  Integer dim, 
  SharedArray<double> densite_initiale, 
  SharedArray<double> energie_initiale, 
  SharedArray<double> pression_initiale, 
  SharedArray<double> temperature_initiale, 
  SharedArray<Real3> vitesse_initiale)  {

}

void EXTERNALService::initUtilisateur(Real3 vitesse_initiale) {  
  Real3 zero = {0.0, 0.0, 0.0};
  pinfo() << vitesse_initiale ;
  if ( vitesse_initiale ==  zero) {
    info() << " Initialisation analytique de " << m_velocity.name();
    ENUMERATE_NODE(inode, allNodes()){
        m_velocity[inode] = - m_node_coord[inode] / m_node_coord[inode].normL2();
        // info() << " noeud " <<  inode.localId() << " " << m_node_coord[inode] << " " << m_velocity[inode];
    }
  } else {
    info() << " Initialisation constante de " << m_velocity.name();
    ENUMERATE_NODE(inode, allNodes()){
        m_velocity[inode] = vitesse_initiale;
    }
  }
}
/*---------------------------------------------------------------------------*/

bool EXTERNALService::hasReverseOption() { return options()->reverseOption;}
Real EXTERNALService::getReverseParameter() { return options()->parametre;}
bool EXTERNALService::isInternalModel() { return false; }

ARCANE_REGISTER_SERVICE_EXTERNAL(EXTERNAL, EXTERNALService); 
