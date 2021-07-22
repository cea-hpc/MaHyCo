// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapALEService.h"
void RemapALEService::ComputeNodeGroupToRelax(){
    
  Int32UniqueArray node_list_lid;
  node_list_lid.clear();
  Integer jj, jp, jm;
  Real3 vecjjp, vecjjm;
  ENUMERATE_NODE(inode, allNodes()) {
    Node node = *inode; 
    Real aircell, sincell, airmax(0.), airmin(1.e10), sinmin(1.);
    ENUMERATE_CELL(icell, node.cells()) {
        Cell cell = * icell;
        for (Integer ii = 0; ii < 4; ++ii) {
            if (cell.node(ii) == node) {
                jj = ii;
                jp = math::abs((jj+1)%4);
                jm = math::abs((jj+3)%4);
            }
        }
        vecjjp = m_node_coord[cell.node(jp)] - m_node_coord[cell.node(jj)];
        vecjjm = m_node_coord[cell.node(jm)] - m_node_coord[cell.node(jj)];
        aircell = math::vecMul2D(vecjjp, vecjjm);
        sincell = aircell / (vecjjp.abs() * vecjjm.abs());
        airmin = ( aircell < airmin ? aircell : airmin );
        airmax = ( aircell > airmax ? aircell : airmax );
        sinmin = ( sincell < sinmin ? sincell : sinmin );
    }
    // pinfo() << node.localId() << " " << airmin/airmax << " " <<  options()->volumCriteria;
    // pinfo() << node.localId() << " " << sinmin << " " <<  options()->angleCriteria;
    if ((airmin/airmax) < options()->volumCriteria) node_list_lid.add(node.localId());
    if (sinmin < options()->angleCriteria) node_list_lid.add(node.localId());
  }
  mesh()->nodeFamily()->createGroup("NodeToRelax", node_list_lid, true);
}
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
  Real3 coordphi, coordpsi, delta;
  Real jacob;
  Real alpha, beta, gamma, weight, dplmax, dplmin;
  Real tauxdlp = 0.00025;
  Real3 cc = {0.5, 0.5, 0.};
  
  
  // sauvegarde de l'ancien maillage 
  m_node_coord_l.copy(m_node_coord);

  NodeGroup Nodes_to_relax = mesh()->nodeFamily()->findGroup("NodeToRelax");
  pinfo() << " nombre de noeuds à relaxer " << Nodes_to_relax.size();
  for( Integer iter=0; iter< options()->nbIterationWinslow ; ++iter){
    ENUMERATE_NODE(inode, Nodes_to_relax){
        Node n1 = *inode;
        if (n1.nbCell() == 4) {
        DirNode dir_nodex(ndmx[inode]);
        Node n6 = dir_nodex.previous();
        Node n2 = dir_nodex.next(); 
        
        DirNode dir_nodey(ndmy[inode]);
        Node n8 = dir_nodey.previous();
        Node n4 = dir_nodey.next();    
        
        DirNode dir_nodexy(ndmy[n2]);
        Node n9 = dir_nodexy.previous();
        Node n3 = dir_nodexy.next();
        
        DirNode dir_nodeyx(ndmy[n6]);
        Node n7 = dir_nodeyx.previous();
        Node n5 = dir_nodeyx.next();
        
        
        coordphi = 0.5*(m_node_coord[n2] - m_node_coord[n6]);
        coordpsi = 0.5*(m_node_coord[n4] - m_node_coord[n8]);
        jacob = coordphi.x * coordpsi.y - coordpsi.x * coordphi.y;
            
                
        alpha = coordphi.abs2();
        beta = 0.5*(coordphi.x * coordpsi.x + coordphi.y * coordpsi.y);
        gamma = coordpsi.abs2();
        weight = 2.*(alpha+gamma);
        if (math::abs(jacob) > 1.e-8 && weight != 0.) {
            delta = (alpha * (m_node_coord[n4] + m_node_coord[n8]) 
                    + gamma * (m_node_coord[n2] + m_node_coord[n6])
            - beta * (m_node_coord[n3] - m_node_coord[n5] + m_node_coord[n7] - m_node_coord[n9]) 
            ) / weight - m_node_coord[n1];
                
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
/**
 *******************************************************************************
 * \file computeVolume
 * \brief Calcul des volumes partiels (lié à chaque face index 0,1,2 et 3 
 *           sur l ancien maillage "m_cell_volume_partial_l"
 *         Calcul des volumes totaux (nouveau "m_cell_volume" 
 *           et ancien maillage "m_cell_volume_l"
 * \param  
 * \return m_cell_volume_partial_l, m_cell_volume_l, m_cell_volume
 *******************************************************************************
 */
void RemapALEService::computeVolumes(){
    
    Real3 u1,u2,u3,u4;
    
    ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        // nouveau volume apres lissage 
        u1 = 0.5 * ( m_node_coord[cell.node(1)] + m_node_coord[cell.node(2)] - m_node_coord[cell.node(0)] - m_node_coord[cell.node(3)]);
        u2 = 0.5 * ( m_node_coord[cell.node(2)] + m_node_coord[cell.node(3)] - m_node_coord[cell.node(0)] - m_node_coord[cell.node(1)]);
        
        m_cell_new_volume[cell] = math::vecMul2D(u1, u2);
        m_cell_new_volume[cell] = m_cell_volume[cell];
        
        // ancien volume issu du lagrange
        u1 = 0.375 * (m_node_coord_l[cell.node(1)]-m_node_coord_l[cell.node(0)]) + 0.125 * (m_node_coord_l[cell.node(2)]-m_node_coord_l[cell.node(3)]);
        u2 = 0.125 * (m_node_coord_l[cell.node(1)]-m_node_coord_l[cell.node(0)]) + 0.375 * (m_node_coord_l[cell.node(2)]-m_node_coord_l[cell.node(3)]);
        
        u3 = 0.375 * (m_node_coord_l[cell.node(3)]-m_node_coord_l[cell.node(0)]) + 0.125 * (m_node_coord_l[cell.node(2)]-m_node_coord_l[cell.node(1)]);
        u4 = 0.125 * (m_node_coord_l[cell.node(3)]-m_node_coord_l[cell.node(0)]) + 0.375 * (m_node_coord_l[cell.node(2)]-m_node_coord_l[cell.node(1)]);
        m_cell_volume_partial_l[cell][0] = math::vecMul2D(u1, u3);
        m_cell_volume_partial_l[cell][1] = math::vecMul2D(u1, u4);
        m_cell_volume_partial_l[cell][2] = math::vecMul2D(u2, u4);
        m_cell_volume_partial_l[cell][3] = math::vecMul2D(u2, u3);
        
        m_cell_volume_l[cell] = m_cell_volume_partial_l[cell][0] + m_cell_volume_partial_l[cell][1] 
                              + m_cell_volume_partial_l[cell][2] + m_cell_volume_partial_l[cell][3];
                              
        // on en profite pour initialiser m_cell_delta_volume
        for (Integer ii = 0 ; ii <  4; ++ii)
            m_cell_delta_volume[cell][ii] = 0.;
    
    }
}
/**
 *******************************************************************************
 * \file computeFlux
 * \brief Calcul des flux aux facex des cellules 
 * \param  
 * \return m_cell_delta_volume[cell][faceindex]
 *******************************************************************************
 */
void RemapALEService::computeFlux(){    
    
    Real deltavol;
    Real3 u1,u2;
    
    FaceDirectionMng fdmx(m_cartesian_mesh->faceDirection(0));
    ENUMERATE_FACE(iface, fdmx.innerFaces()) {
        Face face = *iface;
        DirFace dir_face = fdmx[face];
        Cell cell = dir_face.previousCell();
        Cell cellvois = dir_face.nextCell();
        if (cellvois.localId() == -1 || cell.localId() == -1) continue;
        DirFace dir_facex = fdmx[face];
        u1 = 0.5 * (m_node_coord_l[face.node(1)] + m_node_coord[face.node(1)] - m_node_coord_l[face.node(0)] - m_node_coord[face.node(0)]);
        u2 = 0.5 * (m_node_coord[face.node(0)] + m_node_coord[face.node(1)] - m_node_coord_l[face.node(0)] - m_node_coord_l[face.node(1)]);
        deltavol = math::vecMul2D(u1, u2);
        Integer faceindex(-1), faceindexvois(-1);
        if (dir_facex.previousCell().localId() == cellvois.localId()) {
          faceindex = 3;
          faceindexvois = 1;
        }
        if (dir_facex.nextCell().localId() == cellvois.localId()) {
          faceindex = 1;
          faceindexvois = 3;
        }
        if ( deltavol >= 0.) {
          m_cell_delta_volume[cell][faceindex] = deltavol;
          m_cell_delta_volume[cellvois][faceindexvois] = 0.;
        } else {
          m_cell_delta_volume[cell][faceindex] = 0.;
          m_cell_delta_volume[cellvois][faceindexvois] = -deltavol;
        }
    } 
    FaceDirectionMng fdmy(m_cartesian_mesh->faceDirection(1));
    ENUMERATE_FACE(iface, fdmy.innerFaces()) {
        Face face = *iface;
        DirFace dir_face = fdmy[face];
        Cell cell = dir_face.previousCell();
        Cell cellvois = dir_face.nextCell();
        if (cellvois.localId() == -1 || cell.localId() == -1) continue;
        DirFace dir_facey = fdmy[face];
        u1 = 0.5 * (m_node_coord_l[face.node(1)] + m_node_coord[face.node(1)] - m_node_coord_l[face.node(0)] - m_node_coord[face.node(0)]);
        u2 = 0.5 * (m_node_coord[face.node(0)] + m_node_coord[face.node(1)] - m_node_coord_l[face.node(0)] - m_node_coord_l[face.node(1)]);
        deltavol = math::vecMul2D(u2, u1); // a cause de la numerotation des noeuds des faces 
        Integer faceindex(-1), faceindexvois(-1);
        if (dir_facey.previousCell().localId() == cellvois.localId()) {
          faceindex = 0;
          faceindexvois = 2;
        }
        if (dir_facey.nextCell().localId() == cellvois.localId()) {
          faceindex = 2;
          faceindexvois = 0;
        }
        if ( deltavol >= 0.) {
          m_cell_delta_volume[cell][faceindex] = deltavol;
          m_cell_delta_volume[cellvois][faceindexvois] = 0.;
        } else {
          m_cell_delta_volume[cell][faceindex] = 0.;
          m_cell_delta_volume[cellvois][faceindexvois] = -deltavol;
        }
    }
}
/**
 *******************************************************************************
 * \file computeNewEnvCells
 * \brief Creation des nouvelles envcell et suppressions si besoin
 * \param  
 * \return env_cells
 *******************************************************************************
 */
void RemapALEService::computeNewEnvCells(Integer index_env) {
    
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  Int32UniqueArray cells_to_add;
  Int32UniqueArray cells_to_remove;
  cells_to_add.clear();
  cells_to_remove.clear();
  Integer max_id = allCells().itemFamily()->maxLocalId();
  Int32UniqueArray cells_marker(max_id,-1);
  IMeshEnvironment* env = mm->environments()[index_env];
  ENUMERATE_CELL(icell, env->cells()) {
    cells_marker[icell.localId()] = 0;
  }
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    if (m_phi[cell][index_env] < options()->threshold  && cells_marker[cell.localId()] == 0) {
        cells_to_remove.add(cell.localId());
        pinfo() << " cell " << cell.localId() 
            << " ( " << cell->uniqueId() << " ) " << " fraction " << m_phi[cell][index_env]
            << " retirée dans l'env " << env->name();
    }
        
    if (m_phi[cell][index_env] > options()->threshold  && cells_marker[cell.localId()] == -1) {
        cells_to_add.add(cell.localId());
        pinfo() << " cell " << cell.localId() 
            << " ( " << cell->uniqueId() << " ) "  << " fraction " << m_phi[cell][index_env] 
            << " ajoutée dans l'env apres normalisation " << env->name();
    }
  }
  if (!cells_to_add.empty()) {
    pinfo() << "ADD_CELLS to env " << env->name() << " n=" << cells_to_add.size();
    env->cells().addItems(cells_to_add);
  }
  if (!cells_to_remove.empty()){
    pinfo() << "REMOVE_CELLS to env " << env->name() << " n=" << cells_to_remove.size();
    env->cells().removeItems(cells_to_remove);
  }
}
/**
 *******************************************************************************
 * \file computeApproPhi
 * \brief Calcul approximation d une variable Phi donnees
 * \param  
 * \return m_appro_phi[cell]
 *******************************************************************************
 */
void RemapALEService::computeApproPhi(Integer index_env, 
                                      VariableCellArrayReal CellVolumeOrMassePartiel, 
                                      VariableCellArrayReal DeltaVolumeOrMass){    
    
    Real dx0, dx1, dx2, dx3, dx0p, dx1p, dx2p, dx3p;
    Real dfia0, dfi2a, dfia1, dfi3a;
    Real sign(0), fiprim, fiprimtmp;
    m_appro_phi.fill(0.0);
    CellDirectionMng cdmx(m_cartesian_mesh->cellDirection(0));
     ENUMERATE_CELL(icell, allCells()) {
      Cell cell = * icell;
      DirCell dir_cellx(cdmx.cell(cell));
      Cell cellvois1 = dir_cellx.next();
      Cell cellvois3 = dir_cellx.previous();

      if (DeltaVolumeOrMass[cell][1] == 0. && DeltaVolumeOrMass[cell][3] == 0.) continue;
      fiprim = 0.;
      if (cellvois1.localId() != -1 && cellvois3.localId() != -1) {        
        // deltaX 
        dx1 = CellVolumeOrMassePartiel[cell][1] + 
        CellVolumeOrMassePartiel[cell][2];
        dx3 = CellVolumeOrMassePartiel[cell][3] + 
        CellVolumeOrMassePartiel[cell][0];
        
        dfia1 = m_phi[cellvois1][index_env]  - m_phi[cell][index_env] ; // next - cell
        dfi3a =  m_phi[cell][index_env]  - m_phi[cellvois3][index_env] ; // cell - previous
  
        if ((dfia1*dfi3a) > 0.) {
          sign = 1;
          if ( dfia1 < 0.) {
            sign = -1;
            dfia1 = - dfia1;
            dfi3a = - dfi3a;
          }
          // premiere valeur posible de la derivee fiprim
          fiprim = dfia1/dx1;
          // seconde valeur posible de la derivee fiprim
          fiprimtmp = dfi3a/dx3;
          if ( fiprimtmp < fiprim) fiprim = fiprimtmp;
          dx1p = dx1 + CellVolumeOrMassePartiel[cellvois1][0] + CellVolumeOrMassePartiel[cellvois1][3];
          dx3p = dx3 + CellVolumeOrMassePartiel[cellvois3][2] + CellVolumeOrMassePartiel[cellvois3][1];
          // troisieme valeur posible de la derivee fiprim
          fiprimtmp = (dfi3a*dx1p*dx1p +  dfia1*dx3p*dx3p) / (dx1p*dx3p*(dx1p+dx3p));
          
          if ( fiprimtmp < fiprim) fiprim = fiprimtmp;
          
          fiprim *=sign;
        }
      }
      if (DeltaVolumeOrMass[cell][1] != 0.)
          m_appro_phi[cell][1] = m_phi[cell][index_env]  - fiprim*(dx1-0.5*DeltaVolumeOrMass[cell][1]);
        
      if (DeltaVolumeOrMass[cell][3] != 0.)
          m_appro_phi[cell][3] = m_phi[cell][index_env]  - fiprim*(dx3-0.5*DeltaVolumeOrMass[cell][3]);
    
    }
    CellDirectionMng cdmy(m_cartesian_mesh->cellDirection(1));
    ENUMERATE_CELL(icell, allCells()) {
      Cell cell = * icell;
      DirCell dir_celly(cdmy.cell(cell));
      Cell cellvois0 = dir_celly.previous();
      Cell cellvois2 = dir_celly.next();
      
      if (DeltaVolumeOrMass[cell][0] == 0. && DeltaVolumeOrMass[cell][2] == 0.) continue;
      fiprim = 0.;
      if (cellvois0.localId() != -1 && cellvois2.localId() != -1) {
        // deltaX 
        dx0 = CellVolumeOrMassePartiel[cell][0] + CellVolumeOrMassePartiel[cell][1];
        dx2 = CellVolumeOrMassePartiel[cell][2] + CellVolumeOrMassePartiel[cell][3];
        
        dfia0 = m_phi[cell][index_env]  -  m_phi[cellvois0][index_env] ; // cell - previous 
        dfi2a = m_phi[cellvois2][index_env]  -  m_phi[cell][index_env] ;  // next - cell
        
        
        if ((dfia0*dfi2a) > 0.) {
          sign = 1;
          if ( dfia0 < 0.) {
            sign = -1;
            dfia0 = - dfia0;
            dfi2a = - dfi2a;
          }
          // premiere valeur posible de la derivee fiprim
          fiprim = dfia0/dx0;
          // seconde valeur posible de la derivee fiprim
          fiprimtmp = dfi2a/dx2;
          if ( fiprimtmp < fiprim) fiprim = fiprimtmp;
          // troisieme valeur posible de la derivee fiprim
          dx0p = dx0 + CellVolumeOrMassePartiel[cellvois0][3] + CellVolumeOrMassePartiel[cellvois0][2];
          dx2p = dx2 + CellVolumeOrMassePartiel[cellvois2][0] + CellVolumeOrMassePartiel[cellvois2][1];
          fiprimtmp = (dfi2a*dx0p*dx0p +  dfia0*dx2p*dx2p) / (dx0p*dx2p*(dx0p+dx2p));
          
          if ( fiprimtmp < fiprim) fiprim = fiprimtmp; 
          
          fiprim *=sign;
        }
      }
      if (DeltaVolumeOrMass[cell][0] != 0.)
        m_appro_phi[cell][0] = m_phi[cell][index_env] - fiprim*(dx0-0.5*DeltaVolumeOrMass[cell][0]);
        
      if (DeltaVolumeOrMass[cell][2] != 0.)
        m_appro_phi[cell][2] = m_phi[cell][index_env] - fiprim*(dx2-0.5*DeltaVolumeOrMass[cell][2]);
        
     }   
    
}
/**
 *******************************************************************************
 * \file computeNewPhi
 * \brief Calcul de la nouvelle valeur de Phi dans les nouvelles mailles
 * \param  
 * \return m_phi[cell] 
 *******************************************************************************
 */
void RemapALEService::computeNewPhi(Integer index_env,
                                    VariableCellReal OldVolumeOrMass, 
                                    VariableCellReal NewVolumeOrMass, 
                                    VariableCellArrayReal DeltaVolumeOrMass){   
    
    
    CellDirectionMng cdmx(m_cartesian_mesh->cellDirection(0));
    CellDirectionMng cdmy(m_cartesian_mesh->cellDirection(1));
    ENUMERATE_CELL(icell, allCells()) {
      Cell cell = * icell;
      m_phi[cell][index_env] *= OldVolumeOrMass[cell];
      DirCell dir_cellx(cdmx.cell(cell));
      
      Cell cellvois1 = dir_cellx.next();
      if (cellvois1.localId() != -1)
        m_phi[cell][index_env] += - DeltaVolumeOrMass[cell][1]*m_appro_phi[cell][1] 
                        + DeltaVolumeOrMass[cellvois1][3]*m_appro_phi[cellvois1][3]; 
                        
                        
      Cell cellvois3 = dir_cellx.previous();;
      if (cellvois3.localId() != -1)
        m_phi[cell][index_env] += - DeltaVolumeOrMass[cell][3]*m_appro_phi[cell][3] 
                     + DeltaVolumeOrMass[cellvois3][1]*m_appro_phi[cellvois3][1]; 
                     
      DirCell dir_celly(cdmy.cell(cell));
      
      Cell cellvois0 = dir_celly.previous(); // corection à comprendre
      if (cellvois0.localId() != -1)
        m_phi[cell][index_env] += - DeltaVolumeOrMass[cell][0]*m_appro_phi[cell][0] 
                     + DeltaVolumeOrMass[cellvois0][2]*m_appro_phi[cellvois0][2]; 
      
      Cell cellvois2 = dir_celly.next(); // corection à comprendre
      if (cellvois2.localId() != -1)
        m_phi[cell][index_env] += - DeltaVolumeOrMass[cell][2]*m_appro_phi[cell][2] 
                     + DeltaVolumeOrMass[cellvois2][0]*m_appro_phi[cellvois2][0]; 
      
      m_phi[cell][index_env] /=NewVolumeOrMass[cell];
    } 
}
    
