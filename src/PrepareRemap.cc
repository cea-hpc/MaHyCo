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
 * \return m_face_coord, m_face_length_lagrange,
 *         m_face_normal_velocity
 *******************************************************************************
 */
void MahycoModule::
computeFaceQuantitesForRemap()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans computeFaceQuantitesForRemap()";
  bool csts = options()->schemaCsts();  
  Real one_over_nbnode = m_dimension == 2 ? .5  : .25 ;
  
  if (m_dimension == 3) {
//     ENUMERATE_FACE (iFace, allFaces()) {
//     Face face = *iFace;
//     Real3 face_vec1 = m_node_coord[face.node(2)] - m_node_coord[face.node(0)]; 
//     Real3 face_vec2 = m_node_coord[face.node(3)] - m_node_coord[face.node(1)]; 
//     m_face_length_lagrange[face][0]  = 0.5 * math::abs(produit(face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y));
//     m_face_length_lagrange[face][1]  = 0.5 * math::abs(produit(face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x));
//     m_face_length_lagrange[face][2]  = 0.5 * math::abs(produit(face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x));
//     }
    
    auto fnc = m_connectivity_view.faceNode();
    
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);
    
    auto in_node_coord = ax::viewIn(command,m_node_coord);
    auto out_face_length = ax::viewOut(command,m_face_length_lagrange);
    
    command << RUNCOMMAND_ENUMERATE(Face,fid,allFaces()) {
    
      Real3 coord[4];
      Int32 index=0;
      for( NodeLocalId nid : fnc.nodes(fid) ){
        coord[index]=in_node_coord[nid];
        ++index;
      }
     
      Real3 face_vec1 = coord[2] - coord[0];
      Real3 face_vec2 = coord[3] - coord[1];
      
      Real3 out_face_length_temp = 0.5 * math::cross(face_vec1, face_vec2);
      out_face_length[fid] = out_face_length_temp.absolute();
    };
    
  } 
  else {
     Real3 tempveclen;
     ENUMERATE_FACE (iFace, allFaces()) {
        Face face = *iFace;
        tempveclen = m_node_coord[face.node(1)] - m_node_coord[face.node(0)]; 
        m_face_length_lagrange[face][0] = math::abs(tempveclen.y);
        m_face_length_lagrange[face][1] = math::abs(tempveclen.x);
        m_face_length_lagrange[face][2] = 0.;
    }
  }
  
  if (csts) {
    ENUMERATE_FACE (iFace, allFaces()) {
      Face face = *iFace;
      Real3 vitesse_moy = {0. , 0. , 0.};
      m_face_coord[face]=0.;
      for (Integer inode = 0; inode < face.nbNode(); ++inode) {
        m_face_coord[face] +=  one_over_nbnode * m_node_coord[face.node(inode)];
        vitesse_moy += 0.5 * (m_velocity[face.node(inode)] + m_velocity_n[face.node(inode)]);
      }
      m_face_normal_velocity[face] = math::dot((one_over_nbnode * vitesse_moy), m_face_normal[face]);  
    }
  }
  else {
//     ENUMERATE_FACE (iFace, allFaces()) {
//       Face face = *iFace;
//       Real3 vitesse_moy = {0. , 0. , 0.};
//       m_face_coord[face]=0.;
//       for (Integer inode = 0; inode < face.nbNode(); ++inode) {
//         m_face_coord[face] +=  one_over_nbnode * m_node_coord[face.node(inode)];
//         vitesse_moy += m_velocity[face.node(inode)];
//       }
//       m_face_normal_velocity[face] = math::dot((one_over_nbnode * vitesse_moy), m_face_normal[face]);  
//     }
    
    auto fnc = m_connectivity_view.faceNode();
    
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);
    
    auto in_node_coord = ax::viewIn(command,m_node_coord);
    auto in_velocity = ax::viewIn(command,m_velocity);
    auto in_face_normal = ax::viewIn(command,m_face_normal);
    
    auto out_face_coord = ax::viewOut(command,m_face_coord);
    auto out_face_normal_velocity = ax::viewOut(command,m_face_normal_velocity);
    
    command << RUNCOMMAND_ENUMERATE(Face,fid,allFaces()) {

      Real3 vitesse_moy = {0. , 0. , 0.};
      out_face_coord[fid] = Real3(0., 0., 0.);
      
      Int32 index=0;
      for( NodeLocalId nid : fnc.nodes(fid) ){
        out_face_coord[fid] += one_over_nbnode * in_node_coord[nid];
        vitesse_moy += in_velocity[nid];
        ++index;
      }    
      out_face_normal_velocity[fid] = math::dot((one_over_nbnode * vitesse_moy), in_face_normal[fid]);
    };
  }
  
  m_face_length_lagrange.synchronize();
  m_face_normal_velocity.synchronize();
  PROF_ACC_END;
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
 *                           0 : masse
 *                           1 à 2 : quantite de mouvement
 *                           3 : energie cinetique
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
 *                           0 : densite moyenne
 *                           1 à 2 : vitesse
 *                           3 : energie cinetique specifique
 * \param 
 * \return m_u_lagrange, m_u_dual_lagrange, m_phi_lagrange, m_phi_dual_lagrange
 *******************************************************************************
 */
void MahycoModule::computeVariablesForRemap()
{
  PROF_ACC_BEGIN(__FUNCTION__);
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
      m_u_lagrange[cell][3 * nb_total_env + 0] = 0.;
      m_u_lagrange[cell][3 * nb_total_env + 1] = 0.;
      m_u_lagrange[cell][3 * nb_total_env + 2] = 0.;
      // // energie cinetique centree
      m_u_lagrange[cell][3 * nb_total_env + 3] = 0.;
      // // Pseudo partiel pour la quantité de mouvement
      m_u_lagrange[cell][3 * nb_total_env + 4] = m_cell_volume[ev] * m_pseudo_viscosity[ev];
//       if (cell.localId() == 754) info() << cell.localId() << " pseudo avant proj " << m_pseudo_viscosity[ev] 
//            << " ul " << m_u_lagrange[cell][3 * nb_total_env + 4] << " " << index_env;
      
      if (options()->remap()->hasProjectionPenteBorne() == 1) {     
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
      }
      
    }
  }
  
//   ENUMERATE_NODE(inode, allNodes()){
//     // variables duales
//     // quantité de mouvement
//     m_u_dual_lagrange[inode][0] = m_node_mass[inode] * m_velocity[inode].x;
//     m_u_dual_lagrange[inode][1] = m_node_mass[inode] * m_velocity[inode].y;
//     m_u_dual_lagrange[inode][2] = m_node_mass[inode] * m_velocity[inode].z;
//     // masse nodale    
//     m_u_dual_lagrange[inode][3] = m_node_mass[inode];
//     // projection de l'energie cinétique
//     //     if (options->projectionConservative == 1)
//     m_u_dual_lagrange[inode][4] = 0.5 * m_node_mass[inode] * m_velocity[inode].abs();
// 
//     //         if (limiteurs->projectionAvecPlateauPente == 1) {   
//     // *** variables Phi
//      m_phi_dual_lagrange[inode][0] = m_velocity[inode].x;
//      m_phi_dual_lagrange[inode][1] = m_velocity[inode].y;
//      m_phi_dual_lagrange[inode][2] = m_velocity[inode].z;
//      // masse nodale
//      m_phi_dual_lagrange[inode][3] = m_node_mass[inode];
//      // Phi energie cinétique
//      //     if (options->projectionConservative == 1)
//      m_phi_dual_lagrange[inode][4] = 0.5 * m_velocity[inode].abs();
//   }
  
  auto queue = makeQueue(m_runner);
  auto command = makeCommand(queue);
  
  auto in_node_mass = ax::viewIn(command,m_node_mass);
  auto in_velocity = ax::viewIn(command,m_velocity);
  
  auto out_u_dual_lagrange = ax::viewOut(command,m_u_dual_lagrange);
  auto out_phi_dual_lagrange = ax::viewOut(command,m_phi_dual_lagrange);
  
  command << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
    // variables duales
    // quantité de mouvement
    out_u_dual_lagrange[nid][0] = in_node_mass[nid] * in_velocity[nid].x;
    out_u_dual_lagrange[nid][1] = in_node_mass[nid] * in_velocity[nid].y;
    out_u_dual_lagrange[nid][2] = in_node_mass[nid] * in_velocity[nid].z;
    // masse nodale    
    out_u_dual_lagrange[nid][3] = in_node_mass[nid];
    // projection de l'energie cinétique
    //     if (options->projectionConservative == 1)
    out_u_dual_lagrange[nid][4] = 0.5 * in_node_mass[nid] * in_velocity[nid].normL2();
    
    //         if (limiteurs->projectionAvecPlateauPente == 1) {   
    // *** variables Phi
    out_phi_dual_lagrange[nid][0] = in_velocity[nid].x;
    out_phi_dual_lagrange[nid][1] = in_velocity[nid].y;
    out_phi_dual_lagrange[nid][2] = in_velocity[nid].z;
    // masse nodale
    out_phi_dual_lagrange[nid][3] = in_node_mass[nid];
    // Phi energie cinétique
    //     if (options->projectionConservative == 1)
    out_phi_dual_lagrange[nid][4] = 0.5 * in_velocity[nid].normL2();
  };
  
  PROF_ACC_END;
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
    PROF_ACC_BEGIN(__FUNCTION__);
    debug() << " Entree dans remap()";
    
    Integer withDualProjection = 0;
    // passage des vitesse de n+1/2 à n+1 pour la projection
    if (options()->schemaCsts()) updateVelocityForward();
    if (!options()->sansLagrange) withDualProjection = 1;
    
    computeVariablesForRemap();
    computeFaceQuantitesForRemap();
    
    options()->remap()->appliRemap(m_dimension, withDualProjection, m_nb_vars_to_project, m_nb_env);
    
    // reinitialisaiton des variables (à l'instant N) pour eviter des variables non initialisees
    // pour les nouveaux envirronements crees dans les mailles mixtes par la projection
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
      _computeMultiEnvGlobalCellId(); // la carte des mailles par env a evolué
      computePressionMoyenne();
   }
    PROF_ACC_END;
  }
}
