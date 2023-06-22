#include "CHOCBULLEService.h"

void CHOCBULLEService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void CHOCBULLEService::initVarMono(Integer dim, Real3 densite_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
    
  /* Pas d'implémentation monoMat */
}
void CHOCBULLEService::initMat(Integer dim)  {
    
  
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    double r(0.);
    // bulle surchargera l'aire
    // centre de la bulle
    Real3  Xb = {0.320, 0.089, 0.0};
    //Real3  Xb = {0.040, 0.089, 0.0};
    double rb = 0.025;
    const Real x = m_cell_coord[cell].x - Xb.x;
    const Real y = m_cell_coord[cell].y - Xb.y;
    const Real z = m_cell_coord[cell].z - Xb.z;
    r = std::sqrt(x*x + y*y);
    if (r <= rb) {
      m_materiau[cell] = 0.;
      // m_materiau[cell] = 0.9;
    } else {
      m_materiau[cell] = 1.;
      // m_materiau[cell] = 0.1;
    }
  }
}

void CHOCBULLEService::initVar(Integer dim, Real3 densite_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
 
 // mise à zero puis initialisation des fractions de masses et volumes
 m_mass_fraction.fill(0.0);
 m_fracvol.fill(0.0);
 CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
 ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    double r(0.);
    double pInit;
    double rhoInit;
    // Air partout
    rhoInit = 1.0;
    pInit = 1.e5;
    m_density[cell] = rhoInit;
    m_fracvol[cell] = 1.;
    m_mass_fraction[cell] = 1.;
    m_pressure[cell] = pInit;
    
    // bulle surchargera l'aire
    // centre de la bulle
    Real3  Xb = {0.320, 0.089, 0.0};
    // Real3  Xb = {0.040, 0.089, 0.0};
    double rb = 0.025;
    const Real x = m_cell_coord[cell].x - Xb.x;
    const Real y = m_cell_coord[cell].y - Xb.y;
    const Real z = m_cell_coord[cell].z - Xb.z;
    r = std::sqrt(x*x + y*y);
    
    if (r <= rb) {
      pInit = 1.e5;
      rhoInit = 0.182;
      m_density[cell] = rhoInit;
      m_fracvol[cell] = 1.;
      m_mass_fraction[cell] = 1.;
      m_pressure[cell] = pInit;
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Integer index_env = ev.environmentId();  
            m_density[ev] = rhoInit;
            m_pressure[ev] = pInit;
            // ajout mailles mixtes partouts
            // m_fracvol[ev] = 0.5 * (1-index_env) + 0.5 * index_env;
            // m_mass_fraction[ev] = 0.5 * (1-index_env) + 0.5 * index_env;
        }
      }
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    if (m_node_coord[inode].x >= 0.60) {
        m_velocity[inode] = {-124.824, 0.0, 0.0};
    } else {
        m_velocity[inode] = {0.0, 0.0, 0.0};
    }
  } 
}
/*---------------------------------------------------------------------------*/


bool CHOCBULLEService::hasReverseOption() { return options()->reverseOption;}
Real CHOCBULLEService::getReverseParameter() { return options()->parametre;}

/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_CHOCBULLE(CHOCBULLE, CHOCBULLEService);
