// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"
/**
 *******************************************************************************
 * \file computeFaceQuantitesForRemap()
 * \brief Calcul de quantites aux faces pour la projection :
 *    DxLagrange, du milieu, de la longueur des faces et de leur vitesse normale
 *    avec les coordonées du maillage issues du lagrange
 * 
 * \param  m_cell_coord, m_node_coord, m_face_normal
 *         m_velocity
 * \return m_face_coord, m_face_length_lagrange,
 *         m_face_normal_velocity
 *******************************************************************************
 */
void MahycoModule::
computeFaceQuantitesForRemap()
{
  debug() << " Entree dans computeFaceQuantitesForRemap()";
  
  // calcul su le maillage lagrange de  
  // - la vitesse normale à la face : m_face_normal_velocity
  // - la longueur d'une face : m_face_length_lagrange
  // - les coordonnées de la face de reception : m_face_coord
  
  bool csts = options()->schemaCsts();  
  Real one_over_nbnode = m_dimension == 2 ? .5  : .25 ;
  if (m_dimension == 3) {
    ENUMERATE_FACE (iFace, allFaces()) {
    Face face = *iFace;
    Real3 face_vec1 = m_node_coord[face.node(2)] - m_node_coord[face.node(0)]; 
    Real3 face_vec2 = m_node_coord[face.node(3)] - m_node_coord[face.node(1)]; 
    m_face_length_lagrange[face][0]  = 0.5 * math::abs(produit(face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y));
    m_face_length_lagrange[face][1]  = 0.5 * math::abs(produit(face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x));
    m_face_length_lagrange[face][2]  = 0.5 * math::abs(produit(face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x));
    }
  } else {
     Real3 tempveclen;
     ENUMERATE_FACE (iFace, allFaces()) {
        Face face = *iFace;
        tempveclen = m_node_coord[face.node(1)] - m_node_coord[face.node(0)]; 
        m_face_length_lagrange[face][0] = math::abs(tempveclen.y);
        m_face_length_lagrange[face][1] = math::abs(tempveclen.x);
        m_face_length_lagrange[face][2] = 0.;
    }
  }
  ENUMERATE_FACE (iFace, allFaces()) {
    Face face = *iFace;
    Real3 vitesse_moy = {0. , 0. , 0.};
    m_face_coord[face]=0.;
    for (Integer inode = 0; inode < face.nbNode(); ++inode) {
        m_face_coord[face] +=  one_over_nbnode * m_node_coord[face.node(inode)];
        if (csts)
        vitesse_moy += 0.5 * (m_velocity[face.node(inode)] + m_velocity_n[face.node(inode)]);
        else
        vitesse_moy += m_velocity[face.node(inode)];
        
    }
    m_face_normal_velocity[face] = math::dot((one_over_nbnode * vitesse_moy), m_face_normal[face]);  
  }
  m_face_length_lagrange.synchronize();
  m_face_normal_velocity.synchronize();
  
  if (! options()->remap()->isEuler()) {
    // calcul comme dans InitGeometricValues de 
    // - la normale à la face :  m_face_normal[face] utilisée dans computeGradPhiFace
    // - les coordonnees de la face : m_face_coord 
    // - la normale sortante d'une cellule par une face :  m_outer_face_normal[cell][indexbackface]
    // sur le maillage lissé (m_node_coord_0)
    Real one_over_nbnode = m_dimension == 2 ? .5  : .25 ;
    if ( m_dimension == 3 ) {
        ENUMERATE_FACE ( iFace, allFaces() ) {
            Face face = *iFace;
            Real3 face_vec1 = m_node_coord_0[face.node ( 2 )] - m_node_coord_0[face.node ( 0 )];
            Real3 face_vec2 = m_node_coord_0[face.node ( 3 )] - m_node_coord_0[face.node ( 1 )];
            m_face_normal[iFace].x = produit ( face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y );
            m_face_normal[iFace].y = - produit ( face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x );
            m_face_normal[iFace].z = produit ( face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x );
            m_face_normal[iFace] /= m_face_normal[iFace].normL2();
        }
    } else {
        ENUMERATE_FACE ( iFace, allFaces() ) {
            Face face = *iFace;
            m_face_normal[iFace].x = ( m_node_coord_0[face.node ( 1 )].y - m_node_coord_0[face.node ( 0 )].y );
            m_face_normal[iFace].y = ( m_node_coord_0[face.node ( 1 )].x - m_node_coord_0[face.node ( 0 )].x );
            m_face_normal[iFace] /= m_face_normal[iFace].normL2();
        }
    }
    ENUMERATE_FACE ( iFace, allFaces() ) {
        Face face = *iFace;
        m_face_coord[face] = 0.;
        for ( Integer inode = 0; inode < face.nbNode(); ++inode ) {
            m_face_coord[face] +=  one_over_nbnode * m_node_coord_0[face.node ( inode )];
        }
    }
    ENUMERATE_CELL ( icell,allCells() ) {
        ENUMERATE_FACE ( iface, ( *icell ).faces() ) {
            const Face& face = *iface;
            Integer index = iface.index();
            m_outer_face_normal[icell][index] = ( m_face_coord[face]-m_cell_coord_0[icell] )
                                                / ( m_face_coord[face]-m_cell_coord_0[icell] ).normL2();
        }
    }
    m_face_normal.synchronize();
    m_face_coord.synchronize();
    m_outer_face_normal.synchronize();
  }
}
/**
 *******************************************************************************
 * \file computeVariablesForRemap()
 * \brief Remplissage des variables de la projection et de la projection duale
 *         m_u_lagrange (variables aux mailles)
 *                           de 0 à nbmat-1 : volume partiel,
 *                           de nbmat à 2*nbmat-1 : masse partielle
 *                           de 2*nbmat à 3*nbmat-1 : energie partielle
 *                           de 3*nbmat à 3*nbmat+1 : quantite de mouvement
 *                           3*nbmat+2 : enegie cinetique
 *                           3*nbmat+3 : pseudo-viscosite * volume
 *
 *         m_u_dual_lagrange (variables aux noeuds)
 *                           0 à 2 : quantite de mouvement
 *                           3 : masse nodale
 *                           4 : energie cinetique
 *
 *  Pour l'option projection avec limiteurs pente-borne
 *
 *         m_phi_lagrange (variables aux mailles)
 *                           de 0 à nbmat-1 : fraction volumique
 *                           de nbmat à 2*nbmat-1 : densite partielle
 *                           de 2*nbmat à 3*nbmat-1 : energie specifique
 *partielle de 3*nbmat à 3*nbmat+1 : vitesse 3*nbmat+2 : enegie cinetique
 *specifique 3*nbmat+3 : pseudo-viscosite
 *
 *         m_phi_dual_lagrange (variables aux noeuds)
 *                           0 à 2 : vitesse
 *                           3 : masse nodale
 *                           4 : energie cinetique specifique
 * \param 
 * \return m_u_lagrange, m_u_dual_lagrange, m_phi_lagrange, m_phi_dual_lagrange
 *******************************************************************************
 */
void MahycoModule::computeVariablesForRemap()
{
  debug() << " Entree dans computeVariablesForRemap()";
  Integer nb_total_env = mm->environments().size();
  Integer index_env;

  m_u_lagrange.fill(0.);
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;  
      Cell cell = ev.globalCell();
      index_env = ev.environmentId();
      
      // // volumes matériels (partiels)
      m_u_lagrange[cell][index_env] = m_cell_volume[ev];
      // // masses matériels (partiels)
      m_u_lagrange[cell][nb_total_env + index_env] = m_cell_volume[ev] * m_density[ev];
      // // energies matériels (partiels)
      m_u_lagrange[cell][2 * nb_total_env + index_env] = m_cell_volume[ev] * m_density[ev] * m_internal_energy[ev];
      // // Quantites de mouvement centrees
      //m_u_lagrange[cell][3 * nb_total_env + 0] = 0.;
      //m_u_lagrange[cell][3 * nb_total_env + 1] = 0.;
      //m_u_lagrange[cell][3 * nb_total_env + 2] = 0.;
      // // energie cinetique centree
      //m_u_lagrange[cell][3 * nb_total_env + 3] = 0.;
      // // Pseudo partiel pour la quantité de mouvement
      //m_u_lagrange[cell][3 * nb_total_env + 4] = m_cell_volume[ev] * m_pseudo_viscosity[ev];
      m_u_lagrange[cell][7 * nb_total_env ] = m_cell_volume[cell] * m_pseudo_viscosity[cell];
      /* Indice de phases
      m_u_lagrange[cell][3 * nb_total_env + index_env] = m_cell_volume[ev] * m_frac_phase1[ev];
      m_u_lagrange[cell][4 * nb_total_env + index_env] = m_cell_volume[ev] * m_frac_phase2[ev];
      m_u_lagrange[cell][5 * nb_total_env + index_env] = m_cell_volume[ev] * m_frac_phase3[ev]; */
        
      
      if (options()->remap()->hasProjectionPenteBorne() == 1) {     
        // Cell cell = ev.globalCell();
        m_phi_lagrange[cell][index_env]  = m_fracvol[ev];
        m_phi_lagrange[cell][nb_total_env + index_env] = m_density[ev];
        m_phi_lagrange[cell][2 * nb_total_env + index_env] = m_internal_energy[ev];
        // les phi sur la vitesse et energie cinétique n'existent pas en VnR
        //m_phi_lagrange[cell][3 * nb_total_env + 0] = 0.;
        //m_phi_lagrange[cell][3 * nb_total_env + 1] = 0.;
        //m_phi_lagrange[cell][3 * nb_total_env + 2] = 0.;
        //m_phi_lagrange[cell][3 * nb_total_env + 3] = 0.;
        
        // m_phi_lagrange[cell][3 * nb_total_env + 4] = m_pseudo_viscosity[ev];
        m_phi_lagrange[cell][7 * nb_total_env ] = m_pseudo_viscosity[cell];
        /* Indice de phases
        m_phi_lagrange[cell][3 * nb_total_env + index_env] = m_frac_phase1[ev];
        m_phi_lagrange[cell][4 * nb_total_env + index_env] = m_frac_phase2[ev];
        m_phi_lagrange[cell][5 * nb_total_env + index_env] = m_frac_phase3[ev]; */
      } else {
        for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {
          m_phi_lagrange[cell][ivar] = m_u_lagrange[cell][ivar] / m_cell_volume[cell];
        }
      }
      
    }
  }
  ENUMERATE_NODE(inode, allNodes()){
    // variables duales
    // quantité de mouvement
    m_u_dual_lagrange[inode][0] = m_node_mass[inode] * m_velocity[inode].x;
    m_u_dual_lagrange[inode][1] = m_node_mass[inode] * m_velocity[inode].y;
    m_u_dual_lagrange[inode][2] = m_node_mass[inode] * m_velocity[inode].z;
    // masse nodale    
    m_u_dual_lagrange[inode][3] = m_node_mass[inode];
    // projection de l'energie cinétique
    //     if (options->projectionConservative == 1)
    m_u_dual_lagrange[inode][4] = 0.5 * m_node_mass[inode] * m_velocity[inode].squareNormL2();

    //         if (limiteurs->projectionAvecPlateauPente == 1) {   
    // *** variables Phi
     m_phi_dual_lagrange[inode][0] = m_velocity[inode].x;
     m_phi_dual_lagrange[inode][1] = m_velocity[inode].y;
     m_phi_dual_lagrange[inode][2] = m_velocity[inode].z;
     // masse nodale
     m_phi_dual_lagrange[inode][3] = m_node_mass[inode];
     // Phi energie cinétique
     //     if (options->projectionConservative == 1)
     m_phi_dual_lagrange[inode][4] = 0.5 * m_velocity[inode].squareNormL2();
  }
  
}
/**
 *******************************************************************************
 * \file Remap()
 * \brief Point d'entree Remap
 * appelle pour préparer les variables : 
 *    computeVariablesForRemap() et computeFaceQuantitesForRemap()
 * appelle pour la projection proprement dite
 * 
 **/
void MahycoModule::remap() {
    
  if (options()->withProjection) {
    debug() << " Entree dans remap()";
    
    Integer withDualProjection = 0;
    // passage des vitesse de n+1/2 à n+1 pour la projection
    if (options()->schemaCsts()) updateVelocityForward();
    if (!options()->sansLagrange) withDualProjection = 1;
    
    
    if (! options()->remap()->isEuler()) {
        // lissage du maillage
        info() << "Creation liste des noeuds à relaxer";
        options()->remap()->ComputeNodeGroupToRelax();
        info() << "Deplacement des noeuds : lissage";
        options()->remap()->computeLissage();
        info() << "calcul des volumes du maillage lissé";
        options()->remap()->computeVolumes();
        // TODO à harmoniser 
        // maintenant en ALE
        // - m_node_coord_l contient les coordonnées lagrange - maillage de départ
        // - m_node_coord contient les coordonées du maillage d'arrivée (lissé)
        // en ADI (version euler ou non) 
        // - m_node_coord contient les coordonnées lagrange - maillage de départ
        // - m_node_coord_0 contient les coordonées du maillage d'arrivée (lissé)
        // - m_euler_volume contient le volume des mailles du maillage d'arrivée (lissé)
        
        // pour la projection ADI 
        // calcul des différents termes associées au maillage d'arrivée
        // tels que :
        // - la normale à la face :  m_face_normal[face] utilisée dans computeGradPhiFace
        // - la normale sortante d'une cellule par une face :  m_outer_face_normal[cell][indexbackface]
        // - dirproj (1,0,0) puis (0,1,0) puis (0,0,1) en cartésien 
        // - la vitesse normale à la face : m_face_normal_velocity
        // - la longueur d'une face
        // - les coordonnées de la face de reception : m_face_coord
        
        //  ce sera fait dans computeFaceQuantitesForRemap si options()->getIsEulerScheme() = false
    } 
    
    computeVariablesForRemap();
    computeFaceQuantitesForRemap();
    
    options()->remap()->appliRemap(m_dimension, withDualProjection, m_nb_vars_to_project, m_nb_env);
    
    // reinitialisaiton des variables (à l'instant N) pour eviter des variables non initialisees
    // pour les nouveaux environements crees dans les mailles mixtes par la projection
    m_pseudo_viscosity_n.fill(0.0);
    m_internal_energy_n.fill(0.0);
    m_cell_volume_n.fill(0.0);
    m_pressure_n.fill(0.0);
    m_density_n.fill(0.0);
    m_tau_density.fill(0.0);
    
    m_materiau.fill(0.0);
    for (Integer index_env=0; index_env < m_nb_env ; index_env++) {
        IMeshEnvironment* ienv = mm->environments()[index_env];
        // calcul de la fraction de matiere 
        ENUMERATE_ENVCELL(ienvcell,ienv){
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          m_materiau[cell] += index_env*m_fracvol[ev];
        }
    }
    if (!options()->sansLagrange) {
        for (Integer index_env=0; index_env < m_nb_env ; index_env++) {
        IMeshEnvironment* ienv = mm->environments()[index_env];
        // Calcul de la pression et de la vitesse du son
        options()->environment[index_env].eosModel()->applyEOS(ienv);
      }
      computePressionMoyenne();
   }
  }
}
