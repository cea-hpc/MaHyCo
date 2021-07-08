// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapALEService.h"

Integer RemapALEService::getOrdreProjection() { return options()->ordreProjection;}
bool RemapALEService::hasProjectionPenteBorne() { return options()->projectionPenteBorne;}
bool RemapALEService::hasConservationEnergieTotale() { return options()->conservationEnergieTotale;}
bool RemapALEService::isEuler() {return options()->getIsEulerScheme();}
/**
 *******************************************************************************/
void RemapALEService::appliRemap(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) {
    
    synchronizeUremap();  
    resizeRemapVariables( nb_vars_to_project,  nb_env);
    m_cartesian_mesh = ICartesianMesh::getReference(mesh());
    // m_cartesian_mesh->computeDirections();
    info() << "creation liste des noeuds à relaxer";
    ComputeNodeGroupToRelax();
    
    info() << "deplacement des noeuds : lissage";
    computeLissage();
    ENUMERATE_NODE(inode, allNodes()){
      Node n1 = *inode;/*
      info() << n1.localId() << m_node_coord_l[n1];
      info() << n1.localId() << m_node_coord[n1];*/
    }
    // Calcul des volumes anciens et nouveau et des volumes partiels
    computeVolumes();
    // Calcul des flux de volumes 
    computeFlux();
    
    // en multimateriaux, c'est la qu'il faut creer les nouvelles mailles mixtes
    
    pinfo() << " Projection de la masse " ;
    /************************************************************/
    ENUMERATE_CELL(icell,allCells()){
      Cell c = *icell;
      m_phi[c] = m_density[c];
      m_density_l[c] = m_density[c];
      
    }
    computeApproPhi(m_cell_volume_partial_l, m_cell_delta_volume);
    m_appro_phi.synchronize();
    // calcul de la masse dans les nouvelles cellules
    computeNewPhi(m_cell_volume_l, m_cell_volume, m_cell_delta_volume);
    
    m_node_mass_l.copy(m_node_mass);
    m_node_mass.fill(0.0);
    ENUMERATE_CELL(icell,allCells()){
      Cell c = *icell;
      m_density[c] = m_phi[c];
      m_cell_mass[c] = m_density[c] * m_cell_volume[c];
      Real contrib_node_mass = 0.25 * m_cell_mass[c];
      for( NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode){
        m_node_mass[inode] += contrib_node_mass; 
      }
    }
    if (withDualProjection) {
      // Calcul des masses partiels et des flux de masses utilisant 
      // m_appro_phi qui contient encore la densité approchée aux faces
      ENUMERATE_CELL(icell,allCells()){
       Cell c = *icell;
       m_cell_one[c] = 1.;
       m_cell_zero[c] = 0.;
       for (Integer ii = 0; ii < 4; ++ii) {
        // m_appro_phi contient encore l'approximation de la densité
        m_cell_delta_masse[c][ii] =  m_cell_delta_volume[c][ii] * m_appro_phi[c][ii];
        m_cell_masse_partial_l[c][ii] = m_cell_volume_partial_l[c][ii] * m_density_l[c];
       }
      }
      pinfo() << " Projection de la vitesse en X";
      /************************************************************/
      // Projection de la vitesse X aux mailles
      ENUMERATE_CELL(icell,allCells()){
       Cell c = *icell;
       m_phi[c] = 0.25 * (m_velocity[c.node(0)].x + m_velocity[c.node(1)].x + m_velocity[c.node(2)].x + m_velocity[c.node(3)].x);
      }
      computeApproPhi(m_cell_masse_partial_l, m_cell_delta_masse); // avec masses partiels et des flux de masses
      m_appro_phi.synchronize();
      // mise à zero de la variable phi grace à m_cell_zero
      // calcul de la nouvelle valeur de phi (*m_cell_one, qui vaut 1)
      computeNewPhi(m_cell_zero, m_cell_one, m_cell_delta_masse); 
      
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
       // pinfo() << node.localId() << " AV P " << m_velocity[node].x;
        m_velocity[node].x *= m_node_mass_l[node];
      }
      
      ENUMERATE_CELL(icell, allCells()){  
        Cell c = *icell; 
        for( NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode){
          m_velocity[inode].x += 0.25 * m_phi[icell];
        }
      }
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
        m_velocity[node].x /= m_node_mass[node];
       // pinfo() << node.localId() << " AP P" << m_velocity[node].x;
      }  
      pinfo() << " Projection de la vitesse en Y";
      /************************************************************/
      // Projection de la vitesse Y aux mailles
      ENUMERATE_CELL(icell,allCells()){
       Cell c = *icell;
       m_phi[c] = 0.25 * (m_velocity[c.node(0)].y + m_velocity[c.node(1)].y + m_velocity[c.node(2)].y + m_velocity[c.node(3)].y);
      }
      computeApproPhi(m_cell_masse_partial_l, m_cell_delta_masse); // avec masses partiels et des flux de masses
      m_appro_phi.synchronize();
      // mise à zero de la variable phi grace à m_cell_zero
      // calcul de la nouvelle valeur de phi (*m_cell_one, qui vaut 1)
      computeNewPhi(m_cell_zero, m_cell_one, m_cell_delta_masse); 
      
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
        m_velocity[node].y *= m_node_mass_l[node];
      }
      
      ENUMERATE_CELL(icell, allCells()){  
        Cell c = *icell;
        for( NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode){
          m_velocity[inode].y +=  0.25 * m_phi[icell];
        } 
      }
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
        m_velocity[node].y /= m_node_mass[node];
      }    
       
      /************************************************************/ 
      // Projection de l'energie interne 
      ENUMERATE_CELL(icell,allCells()){
        Cell c = *icell;
        m_phi[c] = m_internal_energy[c] * m_density_l[c];
      }
      computeApproPhi(m_cell_volume_partial_l, m_cell_delta_volume);
      m_appro_phi.synchronize();
      // calcul de la masse dans les nouvelles cellules
      computeNewPhi(m_cell_volume_l, m_cell_volume, m_cell_delta_volume);
        
      ENUMERATE_CELL(icell,allCells()){
        Cell c = *icell;
        m_internal_energy[c] = m_phi[c] / m_density[c];
      }
    
    }
    
}
/**
 *******************************************************************************/
void RemapALEService::resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env) {
    
    
  m_u_lagrange.resize(nb_vars_to_project);
  m_u_dual_lagrange.resize(nb_vars_to_project);
  
  m_phi_dual_lagrange.resize(nb_vars_to_project);
  m_dual_grad_phi.resize(nb_vars_to_project);
  m_grad_phi.resize(nb_vars_to_project);
  m_phi_face.resize(nb_vars_to_project);
  m_grad_phi_face.resize(nb_vars_to_project);
  m_delta_phi_face_ar.resize(nb_vars_to_project);
  m_delta_phi_face_av.resize(nb_vars_to_project);
  m_dual_phi_flux.resize(nb_vars_to_project);
  m_front_flux_mass_env.resize(nb_env);
  m_back_flux_mass_env.resize(nb_env);
  m_cell_volume_partial_l.resize(4);
  m_cell_masse_partial_l.resize(4);
  m_cell_delta_volume.resize(4);
  m_appro_phi.resize(4);
}
/**
 *******************************************************************************/
void RemapALEService::remapVariables(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) {
}
/**
 *******************************************************************************
 * \file synchronizeUremap()
 * \brief phase de synchronisation des variables de projection apres projection
 * \return m_phi_lagrange, m_u_lagrange synchonise sur les mailles fantomes
 *******************************************************************************
 */
void RemapALEService::synchronizeUremap()  {
    debug() << " Entree dans synchronizeUremap()";
    m_density.synchronize();
    m_velocity.synchronize();
    m_internal_energy.synchronize();
    m_est_mixte.synchronize();
    m_est_pure.synchronize();
}/*---------------------------------------------------------------------------*/
ARCANE_REGISTER_SERVICE_REMAPALE(RemapALE, RemapALEService);
/*---------------------------------------------------------------------------*/
