
// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapALEService.h"

/**
 *******************************************************************************
 * \file computeLissage
 * \brief Lissage du maillage avec une methode de winslow 
 *        utilisation du maillage structuré
 * \param  
 * \return m_node_coord
 *******************************************************************************
 */
void RemapALEService::computeLissage(){
    
    
  NodeDirectionMng ndmx(m_cartesian_mesh->nodeDirection(0));
  NodeDirectionMng ndmy(m_cartesian_mesh->nodeDirection(1));
  Real3 coordi, coordj, delta;
  Real jacob;
  Real alpha, beta, gamma, weight, dplmax, dplmin;
  Real tauxdlp = 0.25;
  
  
  // sauvegarde de l'ancien maillage 
  m_node_coord_0.copy(m_node_coord);
  
  for( Integer iter=0; iter< 50; ++iter){
      
   for( Integer i=0; i< mesh()->dimension(); ++i){
    NodeDirectionMng ndm(m_cartesian_mesh->nodeDirection(i));
    ENUMERATE_NODE(inode, ndm.innerNodes() ){
      Node n1 = *inode;
      DirNode dir_nodex(ndmx[inode]);
      Node n9 = dir_nodex.previous();
      Node n5 = dir_nodex.next();
    
      DirNode dir_nodey(ndmy[inode]);
      Node n3 = dir_nodey.previous();
    
      Node n7 = dir_nodey.next();        
      DirNode dir_nodexy(ndmy[n9]);
      Node n2 = dir_nodexy.previous();
      Node n8 = dir_nodexy.next();
    
      DirNode dir_nodeyx(ndmy[n5]);
      Node n4 = dir_nodeyx.previous();
      Node n6 = dir_nodeyx.next();
    
      coordi = 0.5*(m_node_coord[n4] - m_node_coord[n8]);
      coordj = 0.5*(m_node_coord[n2] - m_node_coord[n6]);
      jacob = coordi.x*coordj.y - coordj.x*coordi.y;
        
            
      alpha = coordi.abs2();
      beta = 0.5*(coordi.x * coordj.x + coordi.y * coordj.y);
      gamma = coordj.abs2();
      weight = 2.*(alpha+gamma);
        
      if (math::abs(jacob) > 1.e-8 && weight != 0.) {
        delta = (alpha * (m_node_coord[n2] + m_node_coord[n6]) 
        - beta * (m_node_coord[n3] - m_node_coord[n5] + m_node_coord[n7] - m_node_coord[n9]) 
        + gamma*(m_node_coord[n8] + m_node_coord[n4])) / weight - m_node_coord[n1];
            
        dplmax = tauxdlp*math::min(math::sqrt(alpha),sqrt(gamma));
        dplmin = dplmax/10.;
        
        if ((math::abs(delta.x) > dplmin) || (math::abs(delta.y) > dplmin) || iter ==1) {
            
            delta.x = math::max(delta.x, dplmax);
            delta.x = math::min(delta.x, - dplmax);
            delta.y = math::max(delta.y, dplmax);
            delta.y = math::min(delta.y, - dplmax);
            
            m_node_coord[n1] += delta;
        }
      }
    }
  }
  }
}

void RemapALEService::computeVolumes(){
    
    Real3 u1,u2,u3,u4;
    Real deltavol;
    
    ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        u1 = 0.5 * ( m_node_coord[cell.node(1)] + m_node_coord[cell.node(2)] - m_node_coord[cell.node(0)] - m_node_coord[cell.node(3)]);
        u2 = 0.5 * ( m_node_coord[cell.node(2)] + m_node_coord[cell.node(3)] - m_node_coord[cell.node(0)] - m_node_coord[cell.node(1)]);
        
        m_cell_volume[cell] = math::vecMul2D(u1, u2);
        
        u1 = 0.375 * (m_node_coord_0[cell.node(1)]-m_node_coord_0[cell.node(0)]) + 0.125 * (m_node_coord_0[cell.node(2)]-m_node_coord_0[cell.node(3)]);
        u2 = 0.125 * (m_node_coord_0[cell.node(1)]-m_node_coord_0[cell.node(0)]) + 0.375 * (m_node_coord_0[cell.node(2)]-m_node_coord_0[cell.node(3)]);
        
        u3 = 0.375 * (m_node_coord_0[cell.node(3)]-m_node_coord_0[cell.node(0)]) + 0.125 * (m_node_coord_0[cell.node(2)]-m_node_coord_0[cell.node(1)]);
        u4 = 0.125 * (m_node_coord_0[cell.node(3)]-m_node_coord_0[cell.node(0)]) + 0.375 * (m_node_coord_0[cell.node(2)]-m_node_coord_0[cell.node(1)]);
        m_cell_volume_partial_0[cell][0] = math::vecMul2D(u1, u3);
        m_cell_volume_partial_0[cell][1] = math::vecMul2D(u1, u4);
        m_cell_volume_partial_0[cell][2] = math::vecMul2D(u2, u4);
        m_cell_volume_partial_0[cell][3] = math::vecMul2D(u2, u3);
        
        m_cell_volume_0[cell] = m_cell_volume_partial_0[cell][0] + m_cell_volume_partial_0[cell][1] 
                              + m_cell_volume_partial_0[cell][2] + m_cell_volume_partial_0[cell][3];
    }
    m_cell_delta_volume.fill(0.0);
    ENUMERATE_FACE(iface, allFaces()) {
        const Face& face = *iface;
        Cell cell = face.cell(0);
        Cell cellvois = face.cell(1);
        u1 = 0.5 * (m_node_coord_0[face.node(1)] + m_node_coord[face.node(1)] - m_node_coord_0[face.node(0)] - m_node_coord[face.node(0)]);
        u2 = 0.5 * (m_node_coord[face.node(0)] + m_node_coord[face.node(1)] - m_node_coord_0[face.node(0)] - m_node_coord[face.node(1)]);
        deltavol = math::vecMul2D(u1, u2);
        if ( deltavol >= 0.) {
          m_cell_delta_volume[cell][iface.index()] = deltavol;
          m_cell_delta_volume[cellvois][iface.index()] = 0.;
        } else {
          m_cell_delta_volume[cell][iface.index()] = 0.;
          m_cell_delta_volume[cellvois][iface.index()] = -deltavol;
        }
    }
    
}
