#include "POINTTRIPLEService.h"

void POINTTRIPLEService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void POINTTRIPLEService::initVarMono(Integer dim)  {
    
  // mise à zero puis initialisation des fractions de masses et volumes
  m_mass_fraction.fill(0.0);
  m_fracvol.fill(0.0);
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    double pInit;
    double rhoInit;
    if (m_cell_coord[cell].x <= 0.01) {
      pInit = 1.0;
      rhoInit = 1.0;
    } else {
      if (m_cell_coord[cell].y <= 0.015) {
	pInit = 0.1;
	rhoInit = 1.;
      } else {
	pInit = 0.1;
	rhoInit = 0.1;
      }
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
void POINTTRIPLEService::initMat(Integer dim)  {
    
  info() << options()->casTest;
  if (options()->casTest == TriplePoint) {
        initMatMono(dim);
        return;
  }
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    if (m_cell_coord[cell].x <= 0.01) {
      m_materiau[cell] = 0.;
    } else {
      if (m_cell_coord[cell].y <= 0.015) {
	m_materiau[cell] = 1.;
      } else {
	m_materiau[cell] = 2.;
      }
    }
  }
}

void POINTTRIPLEService::initVar(Integer dim)  {
    
 if (options()->casTest == TriplePoint) {  
   initVarMono(dim);
   return;
 }
 // mise à zero puis initialisation des fractions de masses et volumes
 m_mass_fraction.fill(0.0);
 m_fracvol.fill(0.0);
 CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
 ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    if (m_cell_coord[cell].x <= 0.01) {
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
        }
      }
    } else {
      if (m_cell_coord[cell].y <= 0.015) {
	m_density[cell] = 1.;
	m_pressure[cell] = 0.1;
	m_fracvol[cell] = 1.;
	m_mass_fraction[cell] = 1.;
	if (all_env_cell.nbEnvironment() !=1) {
	  ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
	    EnvCell ev = *ienvcell;  
	    Integer index_env = ev.environmentId();  
            m_density[ev] = 1.;
            m_pressure[ev] = 0.1;
	  }
	}
      } else {
	m_density[cell] = 0.1;
	m_pressure[cell] = 0.1;
	m_fracvol[cell] = 1.;
	m_mass_fraction[cell] = 1.;
	if (all_env_cell.nbEnvironment() !=1) {
	  ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
	    EnvCell ev = *ienvcell;  
	    Integer index_env = ev.environmentId();  
            m_density[ev] = 0.1;
            m_pressure[ev] = 0.1;
	  }
	}
      }
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};
  }
}
/*---------------------------------------------------------------------------*/


bool POINTTRIPLEService::hasReverseOption() { return options()->reverseOption;}
Real POINTTRIPLEService::getReverseParameter() { return options()->parametre;}

/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_POINTTRIPLE(POINTTRIPLE, POINTTRIPLEService);
