#include "INSTARTService.h"

void INSTARTService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void INSTARTService::initVarMono(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
  Real a0 = 0.01;
  Real rhoH = 2.0;
  Real rhoL = 1.0;
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    double pInitH(0.),pInitL(0.);
    double rhoInitH(0.),rhoInitL(0.);
    double X,Y,G,Ymin(10.),Ymax(0.),Xmin(10.),Xmax(0.);
    if (options()->casTest == InstaRTX) {
        X = m_cell_coord[cell].x;
        Y = m_cell_coord[cell].y;
        G = -0.1 ; // options()->gravity.y;
        ENUMERATE_NODE(inode, cell.nodes()) {
            Xmin = std::min(m_node_coord[inode].x, Xmin);
            Xmax = std::max(m_node_coord[inode].x, Xmax);
            Ymin = std::min(m_node_coord[inode].y, Ymin);
            Ymax = std::max(m_node_coord[inode].y, Ymax);
        }
    }
    if (options()->casTest == InstaRTY) {
        X = m_cell_coord[cell].y;
        Y = m_cell_coord[cell].z;
        G = -0.1 ; // options()->gravity.z;
        ENUMERATE_NODE(inode, cell.nodes()) {
            Xmin = std::min(m_node_coord[inode].y, Xmin);
            Xmax = std::max(m_node_coord[inode].y, Xmax);
            Ymin = std::min(m_node_coord[inode].z, Ymin);
            Ymax = std::max(m_node_coord[inode].z, Ymax);
        }
    }
    if (options()->casTest == InstaRTZ) {
        X = m_cell_coord[cell].z;
        Y = m_cell_coord[cell].x;
        G = -0.1 ; // options()->gravity.x;
        ENUMERATE_NODE(inode, cell.nodes()) {
            Xmin = std::min(m_node_coord[inode].z, Xmin);
            Xmax = std::max(m_node_coord[inode].z, Xmax);
            Ymin = std::min(m_node_coord[inode].x, Ymin);
            Ymax = std::max(m_node_coord[inode].x, Ymax);
        }
    }
    Real Ylim = .5 + a0 * cos(6*Pi*X);
    Real fracv(0.);
    // maille coupe Ylim sinon fracb =1 (maille au dessus) ou 0. maille en desous) 
    if ((Ylim < Ymax) && (Ylim > Ymin)) {
        // calcul de la fraction volumique par interpolation linéaire
        // fracv = (Ymax - Ylim) / (Ymax - Ymin);
        // calcul par intégration de la fonction Ylim(X) - Aire jusqu'à la hauteur de la maille
        fracv = (.5 - Ymin) * (Xmax-Xmin) + a0 / (6*Pi) * ( sin(6*Pi*Xmax) - sin(6*Pi*Xmin) );
        fracv /= (Ymax - Ymin)*(Xmax - Xmin); // / volume de la maille   
        fracv = 1. - fracv;  // on cherche la fraction du lourd (au dessus de la ligne Ylim)
    } else if (Ylim >= Ymax) {
        fracv = 0.;  // max de maille en desous
    } else if (Ylim <= Ymin) {
        fracv = 1.;  // min de maille au dessus  
    } 
    // if ( Y > Ylim) {
      pInitH = 1.0 + rhoH * G * (Y-1);
      rhoInitH = rhoH;
    // } else {
      pInitL = 1.0 + rhoH * G *(Ylim-1) 
                  + rhoL * G *(Y-Ylim) ;
      rhoInitL = rhoL;
    // }
    m_density[cell] = rhoInitH*fracv + rhoInitL*(1.-fracv);
    m_pressure[cell] = pInitH*fracv + pInitL*(1.-fracv);
    /* 
     if ((Ylim < Ymax) && (Ylim > Ymin)) {
        pinfo() << "Y " << Y << " et Ylim " << Ylim << " avec Ymin " << Ymin << " et Ymax " << Ymax ;
        pinfo() << " cell " << cell.localId() << " density " << m_density[cell] << " pression " << m_pressure[cell] << " fracv " << fracv;
        }
    */
    
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};
  }
}
void INSTARTService::initMat(Integer dim)  {
    
  info() << options()->casTest;
  if (options()->casTest == InstaRTX ||
       options()->casTest == InstaRTY ||
       options()->casTest == InstaRTZ) {
        initMatMono(dim);
        return;
  }
  Real a0 = 0.01;
  double X,Y;
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    if (options()->casTest == BiInstaRTX) {
        X = m_cell_coord[cell].x;
        Y = m_cell_coord[cell].y;
    }
    if (options()->casTest == BiInstaRTY) {
        X = m_cell_coord[cell].y;
        Y = m_cell_coord[cell].z;
    }
    if (options()->casTest == BiInstaRTZ) {
        X = m_cell_coord[cell].z;
        Y = m_cell_coord[cell].x;
    }
    Real Ylim = .5 + a0 * cos(6*Pi*X); 
    if ( Y > Ylim) {
      m_materiau[cell] = 0.;
    } else {
      m_materiau[cell] = 1.;
    }
  }
}

void INSTARTService::initVar(Integer dim, Real3 densite_initiale, Real3 energie_initiale, Real3 pression_initiale, 
                                   Real3x3 vitesse_initiale)  {
    
    
 if (options()->casTest == InstaRTX ||
       options()->casTest == InstaRTY ||
       options()->casTest == InstaRTZ) {
        initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, vitesse_initiale);
        return;
 }
 
 Real rhoH = 2.0;
 Real rhoL = 1.0;
 Real a0 = 0.01;
 CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
 ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    double pInitH(0.),pInitL(0.);
    double rhoInitH(0.),rhoInitL(0.);
    double X,Y,G,Ymin(10.),Ymax(0.),Xmin(10.),Xmax(0.);
    if (options()->casTest == BiInstaRTX) {
        X = m_cell_coord[cell].x;
        Y = m_cell_coord[cell].y;
        G = -0.1 ; // options()->gravity.y;
        ENUMERATE_NODE(inode, cell.nodes()) {
            Xmin = std::min(m_node_coord[inode].x, Xmin);
            Xmax = std::max(m_node_coord[inode].x, Xmax);
            Ymin = std::min(m_node_coord[inode].y, Ymin);
            Ymax = std::max(m_node_coord[inode].y, Ymax);
        }
    }
    if (options()->casTest == BiInstaRTY) {
        X = m_cell_coord[cell].y;
        Y = m_cell_coord[cell].z;
        G = -0.1 ; // options()->gravity.z;
        ENUMERATE_NODE(inode, cell.nodes()) {
            Xmin = std::min(m_node_coord[inode].y, Xmin);
            Xmax = std::max(m_node_coord[inode].y, Xmax);
            Ymin = std::min(m_node_coord[inode].z, Ymin);
            Ymax = std::max(m_node_coord[inode].z, Ymax);
        }
    }
    if (options()->casTest == BiInstaRTZ) {
        X = m_cell_coord[cell].z;
        Y = m_cell_coord[cell].x;
        G = -0.1 ; // options()->gravity.x;
        ENUMERATE_NODE(inode, cell.nodes()) {
            Xmin = std::min(m_node_coord[inode].z, Xmin);
            Xmax = std::max(m_node_coord[inode].z, Xmax);
            Ymin = std::min(m_node_coord[inode].x, Ymin);
            Ymax = std::max(m_node_coord[inode].x, Ymax);
        }
    }
    Real Ylim = .5 + a0 * cos(6*Pi*X); 
    Real fracv(0.);
    // maille coupe Ylim sinon fracb =1 (maille au dessus) ou 0. maille en desous) 
    if ((Ylim < Ymax) && (Ylim > Ymin)) {
        // calcul de la fraction volumique par interpolation linéaire
        // fracv = (Ymax - Ylim) / (Ymax - Ymin);
        // calcul par intégration de la fonction Ylim(X) - Aire jusqu'à la hauteur de la maille
        fracv = (.5 - Ymin) * (Xmax-Xmin) + a0 / (6*Pi) * ( sin(6*Pi*Xmax) - sin(6*Pi*Xmin) );
        fracv /= (Ymax - Ymin)*(Xmax - Xmin); // / volume de la maille   
        fracv = 1. - fracv;  // on cherche la fraction du lourd (au dessus de la ligne Ylim)
    } else if (Ylim >= Ymax) {
        fracv = 0.;  // max de maille en desous
    } else if (Ylim <= Ymin) {
        fracv = 1.;  // min de maille au dessus  
    }
    // if ( Y > Ylim) {
      pInitH = 1.0 + rhoH * G * (Y-1);
      rhoInitH = rhoH;
    // } else {
      pInitL = 1.0 + rhoH * G *(Ylim-1) 
                  + rhoL * G *(Y-Ylim) ;
      rhoInitL = rhoL;
    // }
    m_density[cell] = rhoInitH*fracv + rhoInitL*(1.-fracv);
    m_pressure[cell] = pInitH*fracv + pInitL*(1.-fracv);
    if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
            EnvCell ev = *ienvcell;  
            if (ev.environmentId() == 0) {
                Integer index_env = ev.environmentId();  
                m_density[ev] = rhoInitH;
                m_pressure[ev] = pInitH;
                m_fracvol[ev] = fracv;
                m_mass_fraction[ev] = m_density[ev] / m_density[cell] * fracv;
            }
            if (ev.environmentId() == 1) {
                Integer index_env = ev.environmentId();  
                m_density[ev] = rhoInitL;
                m_pressure[ev] = pInitL;
                m_fracvol[ev] = 1.-fracv;
                m_mass_fraction[ev] = m_density[ev] / m_density[cell] * (1.-fracv);
            }
        }
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};
  }
}
/*---------------------------------------------------------------------------*/


bool INSTARTService::hasReverseOption() { return options()->reverseOption;}
Real INSTARTService::getReverseParameter() { return options()->parametre;}

/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_INSTART(INSTART, INSTARTService);
