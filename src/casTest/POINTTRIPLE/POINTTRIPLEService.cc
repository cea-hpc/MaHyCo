#include "POINTTRIPLEService.h"

void POINTTRIPLEService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void POINTTRIPLEService::initVarMono(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    double pInit;
    double rhoInit;
    if (m_cell_coord[cell].x <= 0.01) {
      pInit = pression_initiale[0];
      rhoInit = densite_initiale[0];
    } else {
      if (m_cell_coord[cell].y <= 0.015) {
	   pInit = pression_initiale[1];
	   rhoInit = densite_initiale[1];
      } else {
	   pInit = pression_initiale[2];
	   rhoInit = densite_initiale[2];
      }
    }
    m_density[cell] = rhoInit;
    m_pressure[cell] = pInit;
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = vitesse_initiale[0];
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

void POINTTRIPLEService::initVar(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
    
 if (options()->casTest == MonoTriplePoint) {  
   initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, vitesse_initiale);
   return;
 }
 
 Materials::IMeshMaterialMng* mm = IMeshMaterialMng::getReference(mesh());
 CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
 ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    info() << " cell localId "<< cell.localId() << m_cell_coord[cell];
    if (m_cell_coord[cell].x <= 0.01) {
      m_density[cell] = densite_initiale[0];
      m_pressure[cell] = pression_initiale[0];
      info() << " mat 1 cell localId "<< cell.localId() << " " << all_env_cell.nbEnvironment();
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Integer index_env = ev.environmentId();  
          IMeshEnvironment* env = mm->environments()[index_env];
          if (ev.environmentId() == 0) { 
            m_density[ev] = densite_initiale[0];
            m_pressure[ev] = pression_initiale[0];
            info() << " cell localId "<< cell.localId() << " env " << env->name() ;
          }
        }
      }
    } else {
      if (m_cell_coord[cell].y <= 0.015) {
        m_density[cell] = densite_initiale[1];
        m_pressure[cell] = pression_initiale[1];
        info() << " mat 2 cell localId "<< cell.localId() << " " << all_env_cell.nbEnvironment();
        if (all_env_cell.nbEnvironment() !=1) {
            ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                EnvCell ev = *ienvcell;  
                Integer index_env = ev.environmentId(); 
                IMeshEnvironment* env = mm->environments()[index_env];
                if (ev.environmentId() == 1) { 
                    m_density[ev] = densite_initiale[1];
                    m_pressure[ev] = pression_initiale[1];
                }
            }
        }
      } else {
        m_density[cell] = densite_initiale[2];
        m_pressure[cell] = pression_initiale[2];
        info() << " mat 3 cell localId "<< cell.localId() << " " << all_env_cell.nbEnvironment();
        if (all_env_cell.nbEnvironment() !=1) {
            ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                EnvCell ev = *ienvcell;  
                Integer index_env = ev.environmentId();  
                IMeshEnvironment* env = mm->environments()[index_env];
                if (ev.environmentId() == 2) { 
                    m_density[ev] = densite_initiale[2];
                    m_pressure[ev] = pression_initiale[2];
                }
            }
        }
      }
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = vitesse_initiale[0];
  }
  info() << " FIN de initVar";
}
/*---------------------------------------------------------------------------*/


bool POINTTRIPLEService::hasReverseOption() { return options()->reverseOption;}
Real POINTTRIPLEService::getReverseParameter() { return options()->parametre;}
bool POINTTRIPLEService::isInternalModel() { return true; }
void POINTTRIPLEService::initUtilisateur() { }
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_POINTTRIPLE(POINTTRIPLE, POINTTRIPLEService);
