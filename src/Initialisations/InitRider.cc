#include "MahycoModule.h"


void MahycoModule::initMatRiderMono(Real3 Xb)  {
    
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0;
  }
}
void MahycoModule::initMatRider(Real3 Xb)  {
    
  // rayon interne et externe
  double rb(0.15);
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    Real rmin(10.), rmax(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode][0] - Xb[0]) *
                                       (m_node_coord[inode][0] - Xb[0]) +
                                   (m_node_coord[inode][1] - Xb[1]) *
                                       (m_node_coord[inode][1] - Xb[1]));
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
    }
    
    // Air partout
    m_materiau[cell] = 0;
    // bulle surchargera l'aire
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell][0] - Xb[0]) *
                        (m_cell_coord[cell][0] - Xb[0]) +
                    (m_cell_coord[cell][1] - Xb[1]) *
                        (m_cell_coord[cell][1] - Xb[1]));
    
    if (rmax < rb) {
      // maille pure de bulle
    m_materiau[cell] = 1;
    } else if ((rmax >= rb) && (rmin < rb)) {
      double frac_b = (rb - rmin) / (rmax - rmin);
      m_materiau[cell] = (1. - frac_b);
    }
  }
}
void MahycoModule::initVarRiderMono(Real3 Xb)  {
    
  Real3 cc = {0.5, 0.5, 0.};
  // rayon interne et externe
  double rb(0.15);
        
  info() << " boucle sur les mailles";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // parametres maille
    Real rmin(10.), rmax(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode][0] - Xb[0]) *
                                       (m_node_coord[inode][0] - Xb[0]) +
                                   (m_node_coord[inode][1] - Xb[1]) *
                                       (m_node_coord[inode][1] - Xb[1]));
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
    }
    // Air partout
    m_density[cell] = 0.;
    m_pressure[cell] = 0.;
    // bulle surchargera l'aire
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell][0] - Xb[0]) *
                        (m_cell_coord[cell][0] - Xb[0]) +
                    (m_cell_coord[cell][1] - Xb[1]) *
                        (m_cell_coord[cell][1] - Xb[1]));
    if (rmax < rb) {
      // maille pure de bulle
      m_density[cell] = 1.;
      m_pressure[cell] = 0.;
    } else if ((rmax >= rb) && (rmin < rb)) {
      double frac_b = (rb - rmin) / (rmax - rmin);
      m_density[cell] = frac_b;
      m_pressure[cell] = 0.;
    }
  }
  info() << " boucle sur les noeuds";
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};    
    if ( options()->casTest == MonoRiderTx) m_velocity[inode].x = 1.;
    if ( options()->casTest == MonoRiderTy) m_velocity[inode].y = 1.;
    if ( options()->casTest == MonoRiderT45) m_velocity[inode] = {1.0, 1.0, 0.0};  
    if ( options()->casTest == MonoRiderRotation) {
      Real3 dd  = m_node_coord[inode] - cc;
      double theta = std::atan2(dd[1], dd[0]);
      double r = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1]);
      double omega = 4. * Pi;
      m_velocity[inode].x = -r * omega * std::sin(omega * 0. + theta);
      m_velocity[inode].y = r * omega * std::cos(omega * 0. + theta);
    }
    if ( options()->casTest == MonoRiderVortex || 
        options()->casTest == MonoRiderVortexTimeReverse) {
      Real3 dd = m_node_coord[inode];
      m_velocity[inode].x =
          -2. * std::cos(Pi * dd[1]) * std::sin(Pi * dd[1]) *
          std::sin(Pi * dd[0]) * std::sin(Pi * dd[0]);
      m_velocity[inode].y =
          2. * std::cos(Pi * dd[0]) * std::sin(Pi * dd[0]) *
          std::sin(Pi * dd[1]) * std::sin(Pi * dd[1]);
    }
    if ( options()->casTest == MonoRiderDeformation || 
        options()->casTest == MonoRiderDeformationTimeReverse) {
         Real3 dd = m_node_coord[inode] + cc;
          m_velocity[inode].x =
              std::sin(4. * Pi * dd[0]) * std::sin(4. * Pi * dd[1]);
          m_velocity[inode].y  =
              std::cos(4. * Pi * dd[0]) * std::cos(4. * Pi * dd[1]);
        }
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
  info() << " fin de boucle sur les noeuds";
}
void MahycoModule::initVarRider(Real3 Xb)  {
    
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  Real3 cc = {0.5, 0.5, 0.};
  // rayon interne et externe
  double rb(0.15);
        
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // parametres maille
    Real rmin(10.), rmax(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode][0] - Xb[0]) *
                                       (m_node_coord[inode][0] - Xb[0]) +
                                   (m_node_coord[inode][1] - Xb[1]) *
                                       (m_node_coord[inode][1] - Xb[1]));
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
    }
    // Air partout
    m_density[cell] = 0.;
    m_pressure[cell] = 0.;
    // bulle surchargera l'aire
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell][0] - Xb[0]) *
                        (m_cell_coord[cell][0] - Xb[0]) +
                    (m_cell_coord[cell][1] - Xb[1]) *
                        (m_cell_coord[cell][1] - Xb[1]));
    if (rmax < rb) {
      // maille pure de bulle
      m_density[cell] = 1.;
      m_pressure[cell] = 0.;
    } else if ((rmax >= rb) && (rmin < rb)) {
      // cas des cellules mailles mixtes
      double frac_b = (rb - rmin) / (rmax - rmin);
      AllEnvCell all_env_cell = all_env_cell_converter[cell];
      ENUMERATE_CELL_ENVCELL(ienvcell, all_env_cell) {
        EnvCell ev = *ienvcell;
        if (ev.environmentId() == 0) m_density[ev] = 1.;
        if (ev.environmentId() == 1) m_density[ev] = 0.;
        m_pressure[ev] = 0.;
      }
      m_density[cell] = frac_b;
      m_pressure[cell] = 0.;
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};    
    if ( options()->casTest == RiderTx) m_velocity[inode].x = 1.;
    if ( options()->casTest == RiderTy) m_velocity[inode].y = 1.;
    if ( options()->casTest == RiderT45) m_velocity[inode] = {1.0, 1.0, 0.0};  
    if ( options()->casTest == RiderRotation) {
      Real3 dd  = m_node_coord[inode] - cc;
      double theta = std::atan2(dd[1], dd[0]);
      double r = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1]);
      double omega = 4. * Pi;
      m_velocity[inode].x = -r * omega * std::sin(omega * 0. + theta);
      m_velocity[inode].y = r * omega * std::cos(omega * 0. + theta);
    }
    if ( options()->casTest == RiderVortex || 
        options()->casTest == RiderVortexTimeReverse) {
      Real3 dd = m_node_coord[inode];
      m_velocity[inode].x =
          -2. * std::cos(Pi * dd[1]) * std::sin(Pi * dd[1]) *
          std::sin(Pi * dd[0]) * std::sin(Pi * dd[0]);
      m_velocity[inode].y =
          2. * std::cos(Pi * dd[0]) * std::sin(Pi * dd[0]) *
          std::sin(Pi * dd[1]) * std::sin(Pi * dd[1]);
    }
    if ( options()->casTest == RiderDeformation || 
        options()->casTest == RiderDeformationTimeReverse) {
         Real3 dd = m_node_coord[inode] + cc;
          m_velocity[inode].x =
              std::sin(4. * Pi * dd[0]) * std::sin(4. * Pi * dd[1]);
          m_velocity[inode].y  =
              std::cos(4. * Pi * dd[0]) * std::cos(4. * Pi * dd[1]);
    }
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
}

