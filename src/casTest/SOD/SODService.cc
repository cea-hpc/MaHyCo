#include "SODService.h"

void SODService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void SODService::initVarMono(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                    Real3 temperature_initiale, Real3x3 vitesse_initiale)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    double r(0.);
    double pInit;
    double rhoInit;
    double eInit;
    if (options()->casTest == SodCaseX) r = m_cell_coord[cell].x;
    if (options()->casTest == SodCaseY) r = m_cell_coord[cell].y;
    if (options()->casTest == SodCaseZ) r = m_cell_coord[cell].z;
    if (r < 0.5) {
      pInit = pression_initiale[0];
      rhoInit = densite_initiale[0];
    } else {
      pInit = pression_initiale[0]/10.;
      rhoInit = densite_initiale[0]/8.;
    }
    m_density[cell] = rhoInit;
    m_pressure[cell] = pInit;
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = vitesse_initiale[0];
  }
}
void SODService::initMat(Integer dim)  {
    
  info() << options()->casTest;
  if (options()->casTest == SodCaseX ||
       options()->casTest == SodCaseY ||
       options()->casTest == SodCaseZ) {
        initMatMono(dim);
        return;
  }
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    double r(0.);
    if (options()->casTest == BiSodCaseX) r = m_cell_coord[cell].x;
    if (options()->casTest == BiSodCaseY) r = m_cell_coord[cell].y;
    if (options()->casTest == BiSodCaseZ) r = m_cell_coord[cell].z;
    if (options()->casTest == BiSodSph) 
    {
      const Real x = m_cell_coord[cell].x;
      const Real y = m_cell_coord[cell].y;
      const Real z = m_cell_coord[cell].z;
      r = std::sqrt(x*x + y*y + z*z);
    }
    if (r < 0.5) {
      m_materiau[cell] = 0.;
    } else {
      m_materiau[cell] = 1.;
    }
  }
}

void SODService::initVar(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                    Real3 temperature_initiale, Real3x3 vitesse_initiale)  {
    
    
 if (options()->casTest == SodCaseX ||
       options()->casTest == SodCaseY ||
       options()->casTest == SodCaseZ) {
        initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, temperature_initiale, vitesse_initiale);
        return;
 }
 info() << " on rentre ici"; 
 // mise à zero puis initialisation des fractions de masses et volumes
 CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
 ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    double r(0.);
    double eInit;
    if (m_materiau[cell] == 0) {
      m_density[cell] = densite_initiale[0];
      m_pressure[cell] = pression_initiale[0];
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Integer index_env = ev.environmentId();  
            m_density[ev] = densite_initiale[0];
            m_pressure[ev] = pression_initiale[0];
            
        }
      }
    } else {
      m_density[cell] = densite_initiale[1];;
      m_pressure[cell] = pression_initiale[1];
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Integer index_env = ev.environmentId();  
            m_density[ev] = densite_initiale[1];;
            m_pressure[ev] = pression_initiale[1];
        }
      }
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = vitesse_initiale[0];
  }
}
/*---------------------------------------------------------------------------*/


bool SODService::hasReverseOption() { return options()->reverseOption;}
Real SODService::getReverseParameter() { return options()->parametre;}
bool SODService::isInternalModel() { return true; }
void SODService::initUtilisateur() {}

/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SOD(SOD, SODService);
