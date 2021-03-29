// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"
/**
 *******************************************************************************
 * \file computeFaceQuantitesForRemap()
 * \brief Calcul de quantites aux faces pour la projection :
 *    DxLagrange, du milieu, de la longueur des faces et de leur vitesse normale
 *
 * \param  m_cell_coord, m_node_coord, m_face_normal
 *         m_velocity
 * \return m_deltax_lagrange, m_face_coord, m_face_length_lagrange,
 *         m_face_normal_velocity
 *******************************************************************************
 */
void MahycoModule::
computeFaceQuantitesForRemap()
{
  debug() << " Entree dans computeFaceQuantitesForRemap()";
  bool csts = options()->schemaCsts();
  ENUMERATE_FACE (iFace, allCells().innerFaceGroup()) {
    m_deltax_lagrange[iFace] = math::dot(
      (m_cell_coord[iFace->frontCell()] -  m_cell_coord[iFace->backCell()]),
      m_face_normal[iFace]);
  }
  ENUMERATE_FACE (iFace, allFaces()) {
   Face face = *iFace;
   Real3 vitesse_moy = {0. , 0. , 0.};
   Real3 face_vec1 = m_node_coord[face.node(2)] - m_node_coord[face.node(0)]; 
   Real3 face_vec2 = m_node_coord[face.node(3)] - m_node_coord[face.node(1)]; 
   m_face_length_lagrange[face][0]  = 0.5 * math::abs(produit(face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y));
   m_face_length_lagrange[face][1]  = 0.5 * math::abs(produit(face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x));
   m_face_length_lagrange[face][2]  = 0.5 * math::abs(produit(face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x));
   m_face_coord[face]=0.;
   for (Integer inode = 0; inode < 4; ++inode) {
     m_face_coord[face] +=  0.25 * m_node_coord[face.node(inode)];
     if (csts)
       vitesse_moy += 0.5 * (m_velocity[face.node(inode)] + m_velocity_n[face.node(inode)]);
     else
       vitesse_moy += m_velocity[face.node(inode)];
   }
   m_face_normal_velocity[face] = math::dot((0.25 * vitesse_moy), m_face_normal[face]);  
   // info() << "longeur face " << m_face_length_lagrange[face];
   // " vitesse " << m_face_normal_velocity[face] << " et " << m_face_normal[face];
  }
}
/**
 *******************************************************************************
 * \file computeCellQuantitesForRemap()
 * \brief Calcul du centre des mailles pour la projection
 *
 * \param  varlp->XLagrange
 * \return varlp->XcLagrange
 *******************************************************************************
 */
//inline void MahycoModule::
//computeCellQuantitesForRemap()
//{ 

//         varlp->XcLagrange(cCells) =
//             (1.0 / (6.0 * varlp->vLagrange(cCells)) * reduction7);
//       });
//}
/**
 *******************************************************************************
 * \file computeVariablesForRemap()
 * \brief Remplissage des variables de la projection et de la projection duale
 *         varlp->ULagrange (variables aux mailles)
 *                           de 0 à nbmat-1 : volume partiel,
 *                           de nbmat à 2*nbmat-1 : masse partielle
 *                           de 2*nbmat à 3*nbmat-1 : energie partielle
 *                           de 3*nbmat à 3*nbmat+1 : quantite de mouvement
 *                           3*nbmat+2 : enegie cinetique
 *                           3*nbmat+3 : pseudo-viscosite * volume
 *
 *         varlp->UDualLagrange (variables aux noeuds)
 *                           0 : masse
 *                           1 à 2 : quantite de mouvement
 *                           3 : energie cinetique
 *
 *  Pour l'option projection avec limiteurs pente-borne
 *
 *         varlp->Phi (variables aux mailles)
 *                           de 0 à nbmat-1 : fraction volumique
 *                           de nbmat à 2*nbmat-1 : densite partielle
 *                           de 2*nbmat à 3*nbmat-1 : energie specifique
 *partielle de 3*nbmat à 3*nbmat+1 : vitesse 3*nbmat+2 : enegie cinetique
 *specifique 3*nbmat+3 : pseudo-viscosite
 *
 *         varlp->DualPhi (variables aux noeuds)
 *                           0 : densite moyenne
 *                           1 à 2 : vitesse
 *                           3 : energie cinetique specifique
 * \param m_fracvol_env, varlp->vLagrange, m_mass_fraction_env, m_density_nplus1
 *        m_internal_energy_nplus1,
 * \return varlp->ULagrange, varlp->UDualLagrange, varlp->Phi, varlp->DualPhi
 *******************************************************************************
 */
void MahycoModule::
computeVariablesForRemap()
{
  debug() << " Entree dans computeVariablesForRemap()";
  Integer nb_total_env = mm->environments().size();
  Integer index_env;
//   ENUMERATE_CELL(icell, allCells()) {
//     Cell cell = *icell;
//     for (Integer index = 0; index < m_nb_vars_to_project; index++) { 
//       m_u_lagrange[cell][index] = 0.;
//     }
//   }
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
      m_u_lagrange[cell][3 * nb_total_env + 0] = 0.;
      m_u_lagrange[cell][3 * nb_total_env + 1] = 0.;
      m_u_lagrange[cell][3 * nb_total_env + 2] = 0.;
      // // energie cinetique centree
      m_u_lagrange[cell][3 * nb_total_env + 3] = 0.;
      // // Pseudo partiel pour la quantité de mouvement
      m_u_lagrange[cell][3 * nb_total_env + 4] = m_cell_volume[ev] * m_pseudo_viscosity[ev];
      
      if (options()->projectionPenteBorne == 1) {     
        // Cell cell = ev.globalCell();
        m_phi_lagrange[cell][index_env]  = m_fracvol[ev];
        m_phi_lagrange[cell][nb_total_env + index_env] = m_density[ev];
        m_phi_lagrange[cell][2 * nb_total_env + index_env] = m_internal_energy[ev];
        // les phi sur la vitesse et energie cinétique n'existent pas en VnR
        m_phi_lagrange[cell][3 * nb_total_env + 0] = 0.;
        m_phi_lagrange[cell][3 * nb_total_env + 1] = 0.;
        m_phi_lagrange[cell][3 * nb_total_env + 2] = 0.;
        m_phi_lagrange[cell][3 * nb_total_env + 3] = 0.;
        
        m_phi_lagrange[cell][3 * nb_total_env + 4] = m_pseudo_viscosity[ev];
      } else {
        for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {
          m_phi_lagrange[cell][ivar] = m_u_lagrange[cell][ivar] / m_cell_volume[cell];
        }
        //info() << " PREPARE cell " << cell.localId() << " de coord " << m_cell_coord[cell];
        //for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) 
        //  info() << "phi " << ivar <<  "  " << m_phi_lagrange[cell][ivar];
      }
      
      // // m_r_lagrange[ev] = m_density[ev];
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
    m_u_dual_lagrange[inode][4] = 0.5 * m_node_mass[inode] * m_velocity[inode].abs();
//     if ((m_node_coord[inode].x > 0.5 ) && (m_node_coord[inode].x < 0.51 ) && (inode.localId() == 1809 || inode.localId() == 1810) )
//     {
//         info() << inode.localId() << " u_duale prepare " << m_u_dual_lagrange[inode];
//         info() << "vitesse " << m_velocity[inode] << " et masse " << m_node_mass[inode];
//     }
    //         if (limiteurs->projectionAvecPlateauPente == 1) {   
    // *** variables Phi
     m_phi_dual_lagrange[inode][0] = m_velocity[inode].x;
     m_phi_dual_lagrange[inode][1] = m_velocity[inode].y;
     m_phi_dual_lagrange[inode][2] = m_velocity[inode].z;
     // masse nodale
     m_phi_dual_lagrange[inode][3] = m_node_mass[inode];
     // Phi energie cinétique
     //     if (options->projectionConservative == 1)
     m_u_dual_lagrange[inode][4] = 0.5 * m_velocity[inode].abs();
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
    computeVariablesForRemap();
    computeFaceQuantitesForRemap();
    
    IMesh* mesh = defaultMesh();
    Integer nb_dir = mesh->dimension();
    Integer idir(-1);
    m_cartesian_mesh->computeDirections();
    
    for( Integer i=0; i<nb_dir; ++i){
      
      idir = (i + sens_projection)%3;
      // a ameliorer
      String name;
      if (idir == 0) name="FACE_X";
      else if (idir == 1) name="FACE_Y";
      else if (idir == 2) name="FACE_Z";
      debug() << " projection selon la direction " << idir;
      // calcul des gradients des quantites à projeter aux faces 
      computeGradPhiFace(idir, name);
      // calcul des gradients des quantites à projeter aux cellules
      // (avec limiteur ordinaire) 
      // et pour le pente borne, calcul des flux aux faces des cellules
      computeGradPhiCell(idir);
      // calcul de phiFace 
      // qui contient la valeur reconstruite à l'ordre 1, 2 ou 3 des variables projetees 
      // et qui contient les flux des variables projetees avec l'option pente-borne
      computeUpwindFaceQuantitiesForProjection(idir, name);
      
      computeUremap(idir);
      
      if (!options()->sansLagrange) computeDualUremap(idir, name);
    }
    sens_projection++;
    sens_projection = sens_projection%3;
    remapVariables();
  }
}
