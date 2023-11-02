#include "EXTERNALService.h"


void EXTERNALService::initMatMono(Integer dim)  {
}
void EXTERNALService::initMat(Integer dim)  {
   
}
void EXTERNALService::initVarMono(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                    Real3 temperature_initiale, Real3x3 vitesse_initiale)  {
}
void EXTERNALService::initVar(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                    Real3 temperature_initiale, Real3x3 vitesse_initiale)  {

}

void EXTERNALService::initUtilisateur() {  
  info() << " Initialisation de " << m_velocity.name();
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = - m_node_coord[inode] / m_node_coord[inode].normL2();
    // info() << " noeud " <<  inode.localId() << " " << m_node_coord[inode] << " " << m_velocity[inode];
  }
}
/*---------------------------------------------------------------------------*/

bool EXTERNALService::hasReverseOption() { return options()->reverseOption;}
Real EXTERNALService::getReverseParameter() { return options()->parametre;}
bool EXTERNALService::isInternalModel() { return false; }

ARCANE_REGISTER_SERVICE_EXTERNAL(EXTERNAL, EXTERNALService); 
