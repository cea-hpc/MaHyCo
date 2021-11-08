// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapALEService.h"
#include <arcane/ServiceBuilder.h>

/** Constructeur de la classe */
RemapALEService::RemapALEService(const ServiceBuildInfo & sbi)
  : ArcaneRemapALEObject(sbi) {
  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
}
 
Integer RemapALEService::getOrdreProjection() { return options()->ordreProjection;}
bool RemapALEService::hasProjectionPenteBorne() { return options()->projectionPenteBorne;}
bool RemapALEService::hasConservationEnergieTotale() { return options()->conservationEnergieTotale;}
bool RemapALEService::isEuler() {return options()->getIsEulerScheme();}
/**
 *******************************************************************************/
void RemapALEService::appliRemap(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) {
    
    synchronizeUremap();  
    resizeRemapVariables( nb_vars_to_project,  nb_env);
    
    m_cartesian_mesh = CartesianInterface::ICartesianMesh::getReference(mesh());
    mm = IMeshMaterialMng::getReference(mesh());
  
    if (! options()->getIsEulerScheme()) {
      info() << "Creation liste des noeuds à relaxer";
      ComputeNodeGroupToRelax();
      info() << "Deplacement des noeuds : lissage";
      computeLissage();
    } else { 
      // sauvegarde de l'ancien maillage 
      m_node_coord_l.copy(m_node_coord);
      
      pinfo() << " Le lissage consiste à revenir sur le maillage euler";
      // Pour avoir de l'euler 
      m_node_coord.copy(m_node_coord_0);
      //
     Int32UniqueArray node_list_lid;
      ENUMERATE_NODE(inode, allNodes()) {
        Node node= *inode;
        if (node.nbCell() == 4) node_list_lid.add(node.localId());
      }
      mesh()->nodeFamily()->createGroup("NodeToRelax", node_list_lid, true);
    }
    NodeGroup Nodes_to_relax = mesh()->nodeFamily()->findGroup("NodeToRelax");
    
    pinfo() << " Calcul des volumes anciens et nouveau et des volumes partiels";
    // Calcul des volumes anciens et nouveau et des volumes partiels
    computeVolumes();
    pinfo() << " Calcul des flux";
    // Calcul des flux de volumes 
    computeFlux();
    
    
    CellToAllEnvCellConverter all_env_cell_converter(mm);
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      m_fracvol_l[cell] = m_fracvol[cell];  // sauvegarde de la fraction volumique issu du lagrange
      m_density_l[cell] = m_density[cell];
      // m_internal_energy_l[cell] = m_internal_energy[cell];
      // pinfo()  << " AV projection envir " << cell.localId() << " " << m_density_l[cell] << " et " <<  m_fracvol_l[cell];
      AllEnvCell all_env_cell = all_env_cell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) {          
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;      
          m_fracvol_l[ev] = m_fracvol[ev];  // sauvegarde de la fraction volumique issu du lagrange
          m_density_l[ev] = m_density[ev];
        }
      }
    }
    // sauvegarde de la vitesse lagrangienne
    m_velocity_l.copy(m_velocity);
    
    // copie des variables à projeter avant la mise a jour des environnements
     for (Integer index_env=0; index_env < nb_env ; index_env++) { 
        IMeshEnvironment* env = mm->environments()[index_env];
        ENUMERATE_ENVCELL(ienvcell,env){
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          // volume env
          m_phi[cell][index_env] = m_fracvol_l[ev];
          // densité env
          m_phi[cell][nb_env+index_env] = m_density_l[ev] * m_fracvol_l[ev];
          // energie env
          m_phi[cell][2*nb_env+index_env] = m_internal_energy[ev] * m_density_l[ev] * m_fracvol_l[ev];
        }
     }
     ENUMERATE_CELL(icell, allCells()){
       Cell cell = * icell;
       // densite moyenne pour les vitesses
       m_phi[cell][3*nb_env] = m_density_l[cell];
       // Vitesse X
       m_phi[cell][3*nb_env+1] =  0.25 * (m_velocity_l[cell.node(0)].x + m_velocity_l[cell.node(1)].x + 
         m_velocity_l[cell.node(2)].x + m_velocity_l[cell.node(3)].x);
       // Vitesse X
       m_phi[cell][3*nb_env+2] =  0.25 * (m_velocity_l[cell.node(0)].y + m_velocity_l[cell.node(1)].y + 
         m_velocity_l[cell.node(2)].y + m_velocity_l[cell.node(3)].y);
       // energie cinetique
       // m_phi[cell][3*nb_env+3] = 0.5 * m_density_l[cell] * 0.25 * 
       //  (m_velocity_l[cell.node(0)].abs2() + m_velocity_l[cell.node(1)].abs2() 
       // + m_velocity_l[cell.node(2)].abs2() + m_velocity_l[cell.node(3)].abs2());
       m_phi[cell][3*nb_env+3] = 0.5 * m_density_l[cell] * (0.25 * (m_velocity_l[cell.node(0)] + m_velocity_l[cell.node(1)] 
        + m_velocity_l[cell.node(2)] + m_velocity_l[cell.node(3)])).abs2();
     }
    // ajout de la matiere dans les mailles voisines : creation de nouvelles envcell
    if (nb_env >1) { 
     /************************************************************/
     for (Integer index_env=0; index_env < nb_env ; index_env++) { 
        IMeshEnvironment* env = mm->environments()[index_env];
        // pinfo() << " Projection du volume pour l'environement " << env->name() << " nombre de cell " << env->cells().size();

        computeApproPhi(index_env, m_cell_volume_partial_l, m_cell_delta_volume);
        m_appro_phi.synchronize();
        // calcul de la masse dans les nouvelles cellules
        computeNewPhi(index_env, m_cell_volume_l, m_cell_new_volume, m_cell_delta_volume);
        
     
        computeNewEnvCells(index_env);
        
     } 
     // finalisation avant remplissage des variables
     mm->forceRecompute();
     // Ici, la carte des environnements a changé
     m_acc_env->updateMultiEnv(mm);
     
     m_fracvol.fill(0.0);
     
     for (Integer index_env=0; index_env < nb_env ; index_env++) { 
        IMeshEnvironment* env = mm->environments()[index_env];
        ENUMERATE_ENVCELL(ienvcell,env){
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          m_fracvol[ev] = m_phi[cell][index_env];
          AllEnvCell all_env_cell = all_env_cell_converter[cell];
          if (all_env_cell.nbEnvironment() !=1) m_fracvol[cell] += m_fracvol[ev];
        }
     }
     // nomalisation
     ENUMERATE_ENV(ienv,mm){
        IMeshEnvironment* env = *ienv;
        ENUMERATE_ENVCELL(ienvcell,env){
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          m_fracvol[ev] /= m_fracvol[cell];
          m_cell_volume[ev] = m_fracvol[ev] * m_cell_volume[cell];
        }
     }
    }  
    // il reste à mettre fravol[cell] à 1 et à remettre m_phi à zero
    ENUMERATE_CELL(icell,allCells()){
      Cell c = *icell;
      for (Integer index_env=0; index_env < nb_env ; index_env++) 
        m_phi[c][index_env] = 0.;
      m_fracvol[c] =1.;  
    }
    // 
    
    m_appro_phi.fill(0.0);
    pinfo() << " Projection de la masse par envirronement " ;
    /************************************************************/
    for (Integer index_env=0; index_env < nb_env ; index_env++) { 
      IMeshEnvironment* env = mm->environments()[index_env];
      pinfo() << " Projection de la masse pour l'envirronement " << env->name();
    
      
      computeApproPhi(nb_env+index_env, m_cell_volume_partial_l, m_cell_delta_volume);
      m_appro_phi.synchronize();
      // calcul de la masse dans les nouvelles cellules
      computeNewPhi(nb_env+index_env, m_cell_volume_l, m_cell_new_volume, m_cell_delta_volume);
    
      // pour la projection de la quantite de mouvement : approxiation de la denstité moyenne
      m_appro_density.fill(0.0);
      ENUMERATE_CELL(icell,allCells()) {
        Cell c = *icell;
        for (Integer ii = 0; ii < 4; ++ii) 
          m_appro_density[c][ii] = m_appro_phi[c][ii];
      }
      
      ENUMERATE_ENVCELL(ienvcell,env){
        EnvCell ev = *ienvcell;
        Cell cell = ev.globalCell();
        m_density[ev] = m_phi[cell][nb_env+index_env] / m_fracvol[ev];
        // pinfo() << env->name() << " projection densité " << cell.localId() << " " << m_density[ev] 
        // << " et " <<  m_fracvol[ev] << " appro phi " << m_appro_density[cell];
      }
      
      m_appro_phi.fill(0.0);
      pinfo() << " Projection de l'energie pour l'envirronement " << env->name();
      
      computeApproPhi(2*nb_env+index_env, m_cell_volume_partial_l, m_cell_delta_volume);
      m_appro_phi.synchronize();
      // calcul des energie massique dans les nouvelles cellules
      computeNewPhi(2*nb_env+index_env, m_cell_volume_l, m_cell_new_volume, m_cell_delta_volume);
        
      ENUMERATE_ENVCELL(ienvcell,env){
        EnvCell ev = *ienvcell;
        Cell cell = ev.globalCell();
        m_internal_energy[ev] = m_phi[cell][2*nb_env+index_env] / m_phi[cell][nb_env+index_env];  // / m_density[ev] / m_fracvol[ev];
      }
    }
    m_appro_phi.fill(0.0);
    if (nb_env >1) { 
      CellToAllEnvCellConverter all_env_cell_converter(mm);
      ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if (all_env_cell.nbEnvironment() !=1) {
          m_density[cell] = 0.;
          ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
            EnvCell ev = *ienvcell;      
            m_density[cell] += m_fracvol[ev] * m_density[ev];
          }
        }
      }
      // calcul de la fraction massique alpha_i*rho_i = c_i * rho
      // alpha_i*rho_i est la quantite projete : m_phi[cell][nb_env+index_env]
      for (Integer index_env=0; index_env < nb_env ; index_env++) { 
        IMeshEnvironment* env = mm->environments()[index_env];
        ENUMERATE_ENVCELL(ienvcell,env) {
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          m_mass_fraction[ev] = m_phi[cell][nb_env+index_env] / m_density[cell];
        }
      }
      ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if (all_env_cell.nbEnvironment() !=1) {
          m_internal_energy[cell] = 0.;
          ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
            EnvCell ev = *ienvcell;      
            m_internal_energy[cell] += m_mass_fraction[ev] * m_internal_energy[ev];
          }
        }
      }
    }
    m_node_mass_l.copy(m_node_mass);
    m_node_mass.fill(0.0);
    ENUMERATE_CELL(icell,allCells()){
      Cell c = *icell;
      AllEnvCell all_env_cell = all_env_cell_converter[c];
      m_cell_mass[c] = m_density[c] * m_cell_volume[c];
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;      
          m_cell_mass[ev] = m_mass_fraction[ev] * m_cell_mass[c];
        }
      }
      Real contrib_node_mass = 0.25 * m_cell_mass[c];
      for( NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode){
        m_node_mass[inode] += contrib_node_mass; 
      }
    }

    if (withDualProjection) {
      if (nb_env >1) { 
        m_appro_phi.fill(0.0);
        pinfo() << " Calcul de l'approximation de la densité aux faces m_appro_density " 
                << " pour la projection de la quantité de mouvement" ;
        // on refait la projection de la densité totale pour calculer appro_phi 
        /************************************************************/
 
        computeApproPhi(3*nb_env, m_cell_volume_partial_l, m_cell_delta_volume);
        m_appro_phi.synchronize();
        // calcul des energie massique dans les nouvelles cellules
        computeNewPhi(3*nb_env, m_cell_volume_l, m_cell_new_volume, m_cell_delta_volume);
        m_appro_density.fill(0.0);
        ENUMERATE_CELL(icell,allCells()) {
            Cell c = *icell;
            for (Integer ii = 0; ii < 4; ++ii) 
            m_appro_density[c][ii] = m_appro_phi[c][ii];
        }
      }
      // Calcul des masses partiels et des flux de masses utilisant 
      // m_appro_density qui contient la densité approchée aux faces
      ENUMERATE_CELL(icell,allCells()){
       Cell c = *icell;
       m_cell_one[c] = 1.;
       m_cell_zero[c] = 0.;
       for (Integer ii = 0; ii < 4; ++ii) {
        // m_appro_density contient  l'approximation de la densité
        m_cell_delta_masse[c][ii] =  m_cell_delta_volume[c][ii] * m_appro_density[c][ii];
        m_cell_masse_partial_l[c][ii] = m_cell_volume_partial_l[c][ii] * m_density_l[c];
       }
      }
      
      m_appro_phi.fill(0.0);
      pinfo() << " Projection de la vitesse en X";
      /************************************************************/
      // Projection de la vitesse X aux mailles
 
      computeApproPhi(3*nb_env+1, m_cell_masse_partial_l, m_cell_delta_masse); // avec masses partiels et des flux de masses
      m_appro_phi.synchronize();
      // mise à zero de la variable phi grace à m_cell_zero
      // calcul de la nouvelle valeur de phi (*m_cell_one, qui vaut 1)
      computeNewPhi(3*nb_env+1, m_cell_zero, m_cell_one, m_cell_delta_masse); 
      
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
        m_velocity[node].x *= m_node_mass_l[node];
      }
      
      ENUMERATE_CELL(icell, allCells()){  
        Cell c = *icell; 
        for( NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode){
          m_velocity[inode].x += 0.25 * m_phi[icell][3*nb_env+1];
        }
      }
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
        m_velocity[node].x /= m_node_mass[node];
      }  
      m_appro_phi.fill(0.0);
      pinfo() << " Projection de la vitesse en Y";
      /************************************************************/
      // Projection de la vitesse Y aux mailles
 
      computeApproPhi(3*nb_env+2, m_cell_masse_partial_l, m_cell_delta_masse); // avec masses partiels et des flux de masses
      m_appro_phi.synchronize();
      // On peut utiliser computeNewPhi avec
      // la mise à zero de la variable phi grace à m_cell_zero
      // et le calcul de la nouvelle valeur de phi est ensuite *m_cell_one (qui vaut 1) donc ne change pas sa valeur.
      computeNewPhi(3*nb_env+2, m_cell_zero, m_cell_one, m_cell_delta_masse); 
      
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
        m_velocity[node].y *= m_node_mass_l[node];
      }
      
      ENUMERATE_CELL(icell, allCells()){  
        Cell c = *icell;
        for( NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode){
          m_velocity[inode].y +=  0.25 * m_phi[icell][3*nb_env+2];
        } 
      }
      ENUMERATE_NODE(inode, allNodes()){
        Node node= *inode;
        m_velocity[node].y /= m_node_mass[node];
      }    
      pinfo() << " Fin de la Projection";
      
      // récuperation du delta d'nergie cinétique en energie interne
      if (hasConservationEnergieTotale()) {
        /************************************************************/ 
        pinfo() << " Projection de l'energie cinétique";
 
        computeApproPhi(3*nb_env+3, m_cell_volume_partial_l, m_cell_delta_volume);
        m_appro_phi.synchronize();
        // calcul de la masse dans les nouvelles cellules
        computeNewPhi(3*nb_env+3, m_cell_volume_l, m_cell_new_volume, m_cell_delta_volume);
        
        ENUMERATE_CELL(icell,allCells()){
          Cell c = *icell;
          Real ec_proj(0.);
          Real delta_ec(0.);
          Real ec_reconst(0.);
          
          ec_proj = m_phi[c][3*nb_env+3] / m_density[c];
          // ec_reconst = 0.5 * 0.25 * 
          //   (m_velocity[c.node(0)].abs2() + m_velocity[c.node(1)].abs2() + m_velocity[c.node(2)].abs2() + m_velocity[c.node(3)].abs2());
          ec_reconst = 0.5 * (0.25 *(m_velocity[c.node(0)] + m_velocity[c.node(1)] + m_velocity[c.node(2)] + m_velocity[c.node(3)])).abs2();
          delta_ec = std::max( 0., ec_proj - ec_reconst);
          m_internal_energy[c] += delta_ec;
        }
    
      }
    }
}
/**
 *******************************************************************************/
void RemapALEService::resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env) {

  m_cell_volume_partial_l.resize(4);
  m_cell_masse_partial_l.resize(4);
  m_cell_delta_volume.resize(4);
  m_cell_delta_masse.resize(4);
  m_appro_phi.resize(4);
  m_appro_density.resize(4);
  m_phi.resize(nb_vars_to_project);
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
