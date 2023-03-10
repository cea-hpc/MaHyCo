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
  if (options()->casTest == MonoTriplePoint) {
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
    
 if (options()->casTest == MonoTriplePoint) {  
   initVarMono(dim);
   return;
 }
 
 Materials::IMeshMaterialMng* mm = IMeshMaterialMng::getReference(mesh());
 info() << " Ici : Cas test " << options()->casTest << " allCells() " << allCells().size();
 // mise à zero puis initialisation des fractions de masses et volumes
 m_mass_fraction.fill(0.0);
 m_fracvol.fill(0.0);
 CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
 ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    info() << " cell localId "<< cell.localId() << m_cell_coord[cell];
    if (m_cell_coord[cell].x <= 0.01) {
      m_density[cell] = 1.0;
      m_pressure[cell] = 1.0;
      m_fracvol[cell] = 1.;
      m_mass_fraction[cell] = 1.;
      info() << " mat 1 cell localId "<< cell.localId() << " " << all_env_cell.nbEnvironment();
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Integer index_env = ev.environmentId();  
          IMeshEnvironment* env = mm->environments()[index_env];
          if (ev.environmentId() == 0) { 
            m_density[ev] = 1.0;
            m_pressure[ev] = 1.0;
            info() << " cell localId "<< cell.localId() << " env " << env->name() ;
          }
        }
      }
    } else {
      if (m_cell_coord[cell].y <= 0.015) {
        m_density[cell] = 1.;
        m_pressure[cell] = 0.1;
        m_fracvol[cell] = 1.;
        m_mass_fraction[cell] = 1.;
        info() << " mat 2 cell localId "<< cell.localId() << " " << all_env_cell.nbEnvironment();
        if (all_env_cell.nbEnvironment() !=1) {
            ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                EnvCell ev = *ienvcell;  
                Integer index_env = ev.environmentId(); 
                IMeshEnvironment* env = mm->environments()[index_env];
                if (ev.environmentId() == 1) { 
                    m_density[ev] = 1.;
                    m_pressure[ev] = 0.1;
                    info() << " cell localId "<< cell.localId() << " env " << env->name() ;
                }
            }
        }
      } else {
        m_density[cell] = 0.1;
        m_pressure[cell] = 0.1;
        m_fracvol[cell] = 1.;
        m_mass_fraction[cell] = 1.;
        info() << " mat 3 cell localId "<< cell.localId() << " " << all_env_cell.nbEnvironment();
        if (all_env_cell.nbEnvironment() !=1) {
            ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                EnvCell ev = *ienvcell;  
                Integer index_env = ev.environmentId();  
                IMeshEnvironment* env = mm->environments()[index_env];
                if (ev.environmentId() == 2) { 
                    m_density[ev] = 0.1;
                    m_pressure[ev] = 0.1;
                    info() << " cell localId "<< cell.localId() << " env " << env->name() ;
                }
            }
        }
      }
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};
  }
  info() << " FIN de initVar";
}
/*---------------------------------------------------------------------------*/


bool POINTTRIPLEService::hasReverseOption() { return options()->reverseOption;}
Real POINTTRIPLEService::getReverseParameter() { return options()->parametre;}

/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_POINTTRIPLE(POINTTRIPLE, POINTTRIPLEService);
