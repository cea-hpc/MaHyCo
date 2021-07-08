#include "SODService.h"

void SODService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void SODService::initVarMono(Integer dim)  {
    
  // mise à zero puis initialisation des fractions de masses et volumes
  m_mass_fraction.fill(0.0);
  m_fracvol.fill(0.0);
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
      pInit = 1.0;
      rhoInit = 1.0;
    } else {
      pInit = 0.1;
      rhoInit = 0.125;
    }
    m_density[cell] = rhoInit;
    m_pressure[cell] = pInit;
    m_fracvol[cell] = 1.;
    m_mass_fraction[cell] = 1.;
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};
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
    if (r < 0.5) {
      m_materiau[cell] = 0.;
      // m_materiau[cell] = 0.9;
    } else {
      m_materiau[cell] = 1.;
      // m_materiau[cell] = 0.1;
    }
  }
}

void SODService::initVar(Integer dim)  {
    
    
 if (options()->casTest == SodCaseX ||
       options()->casTest == SodCaseY ||
       options()->casTest == SodCaseZ) {
        initVarMono(dim);
        return;
 }
 pinfo() << " on rentre ici"; 
 // mise à zero puis initialisation des fractions de masses et volumes
 m_mass_fraction.fill(0.0);
 m_fracvol.fill(0.0);
 CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
 ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    double r(0.);
    double eInit;
    if (options()->casTest == BiSodCaseX) r = m_cell_coord[cell].x;
    if (options()->casTest == BiSodCaseY) r = m_cell_coord[cell].y;
    if (options()->casTest == BiSodCaseZ) r = m_cell_coord[cell].z;
    if (r < 0.5) {
      m_density[cell] = 1.0;
      m_pressure[cell] = 1.0;
      m_fracvol[cell] = 1.;
      m_mass_fraction[cell] = 1.;
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Integer index_env = ev.environmentId();  
            m_density[ev] = 1.0;
            m_pressure[ev] = 1.0;
            // ajout mailles mixtes partouts
            // m_fracvol[ev] = 0.5 * (1-index_env) + 0.5 * index_env;
            // m_mass_fraction[ev] = 0.5 * (1-index_env) + 0.5 * index_env;
            
        }
      }
    } else {
      m_density[cell] = 0.125;
      m_pressure[cell] = 0.1;
      m_fracvol[cell] = 1.;
      m_mass_fraction[cell] = 1.;
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Integer index_env = ev.environmentId();  
            m_density[ev] = 0.125;
            m_pressure[ev] = 0.1;
            // ajout mailles mixtes partouts
            // m_fracvol[ev] = 0.5 * index_env + 0.5 * (1-index_env);
            // m_mass_fraction[ev] = 0.5 * index_env + 0.5 * (1-index_env);
        }
      }
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};
  }
}
/*---------------------------------------------------------------------------*/


bool SODService::hasReverseOption() { return options()->reverseOption;}
Real SODService::getReverseParameter() { return options()->parametre;}

/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SOD(SOD, SODService);
