// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapALEService.h"

Integer RemapALEService::getOrdreProjection() { return options()->ordreProjection;}
bool RemapALEService::hasProjectionPenteBorne() { return options()->projectionPenteBorne;}
bool RemapALEService::hasConservationEnergieTotale() { return options()->conservationEnergieTotale;}
/**
 *******************************************************************************/
void RemapALEService::appliRemap(Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) {
    
    synchronizeUremap();  
    Integer idir = 0; // pas de direction, ce n'est pas de l'ADI
    // deplacement des noeuds : lissage
    computeLissage();
      
    // calcul des gradients des quantites à projeter aux faces 
    computeGradPhiFace(idir, nb_vars_to_project, nb_env);
    // calcul des gradients des quantites à projeter aux cellules
    // (avec limiteur ordinaire) 
    // et pour le pente borne, calcul des flux aux faces des cellules
    computeGradPhiCell(idir, nb_vars_to_project, nb_env);
    // calcul de m_phi_face
    // qui contient la valeur reconstruite à l'ordre 1, 2 ou 3 des variables projetees 
    // et qui contient les flux des variables projetees avec l'option pente-borne
    computeUpwindFaceQuantitiesForProjection(idir, nb_vars_to_project, nb_env);
    
    
    computeUremap(idir, nb_vars_to_project, nb_env);
    synchronizeUremap();
      
//       if (withDualProjection) {
//         computeDualUremap(idir, nb_env);
//         synchronizeDualUremap();
//       }
    
}
/**
 *******************************************************************************/
void RemapALEService::resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env) {
    
    
  m_u_lagrange.resize(nb_vars_to_project);
  m_u_dual_lagrange.resize(nb_vars_to_project);
  m_phi_lagrange.resize(nb_vars_to_project);
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
  m_cell_volume_partial.resize(4);
  m_cell_delta_volume.resize(4);
}
/**
 *******************************************************************************
 * \file computeGradPhiFace1()
 * \brief phase de projection : premiere etapes
 *  calcul des gradients aux faces dans le sens de la projection
 *  calcul des largeurs de cellule dans le sens de la projection
 * \param \return m_grad_phi, m_h_cell_lagrange
 *******************************************************************************
 */
void RemapALEService::computeGradPhiFace(Integer idir, Integer nb_vars_to_project, Integer nb_env)  {
  debug() << " Entree dans computeGradPhiFace()";
  m_h_cell_lagrange.fill(0.0);
  
  FaceDirectionMng fdm(m_cartesian_mesh->faceDirection(idir));
  ENUMERATE_FACE(iface, fdm.allFaces()) {
    Face face = *iface; 
    m_is_dir_face[face][idir] = true;
  }
  if (options()->ordreProjection > 1) {
    ENUMERATE_FACE(iface, fdm.innerFaces()) {
      Face face = *iface; 
      DirFace dir_face = fdm[face];
      Cell cellb = dir_face.previousCell();
      Cell cellf = dir_face.nextCell();
      m_deltax_lagrange[face] = math::dot(
       (m_cell_coord[cellf] -  m_cell_coord[cellb]), m_face_normal[face]); 

      for (Integer ivar = 0 ; ivar <  nb_vars_to_project ; ++ivar) {
        m_grad_phi_face[iface][ivar] = (m_phi_lagrange[cellf][ivar] - m_phi_lagrange[cellb][ivar]) 
                                    / m_deltax_lagrange[iface];
      }
      // somme des distances entre le milieu de la maille et le milieu de la face
      m_h_cell_lagrange[cellb] +=  (m_face_coord[iface] - m_cell_coord[cellb]).abs();
      m_h_cell_lagrange[cellf] +=  (m_face_coord[iface] - m_cell_coord[cellf]).abs();     
    }
  }
  m_grad_phi_face.synchronize();
  m_h_cell_lagrange.synchronize();
}
/**
 *******************************************************************************
 * \file computeGradPhi()
 * \brief phase de projection : seconde etape
 *        calcul du gradient aux mailles limites : m_grad_phi_face
 *        calcul des flux pente-borne arriere ou avant aux mailles
 *        à partir du gradient precedent : m_delta_phi_face_ar, m_delta_phi_face_av
 * \param
 * \return m_grad_phi_face, m_delta_phi_face_ar, m_delta_phi_face_av
 *******************************************************************************
 */
void RemapALEService::computeGradPhiCell(Integer idir, Integer nb_vars_to_project, Integer nb_env) {
    
 
}
/**
*****************************************************************************
 * \file computeUpwindFaceQuantitiesForProjection()
 * \brief  phase de projection : troisieme etap
 *        calcul de m_phi_face
 *     qui contient la valeur reconstruite à l'ordre 1, 2 ou 3 des variables
 * projetees qui contient les flux des variables projetees avec l'option
 * pente-borne \param 
 * \return m_phi_face
*******************************************************************************
*/
void RemapALEService::computeUpwindFaceQuantitiesForProjection(Integer idir, Integer nb_vars_to_project, Integer nb_env) {
    
 
}
/**
 *******************************************************************************
 * \file computeUremap()
 * \brief phase de projection : etape finale
 *        calcul de la variable Uremap1 par ajout ou retrait des flux
          Mises à jour de l'indicateur mailles mixtes
          calcul de la valeur de Phi issu de Uremap1
 * \param
 * \return m_u_lagrange, m_phi_lagrange, m_est_mixte, m_est_pure
 *******************************************************************************
 */
void RemapALEService::computeUremap(Integer idir, Integer nb_vars_to_project, Integer nb_env)  {
    

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
    m_phi_lagrange.synchronize();
    m_u_lagrange.synchronize();
    m_est_mixte.synchronize();
    m_est_pure.synchronize();
}/*---------------------------------------------------------------------------*/
ARCANE_REGISTER_SERVICE_REMAPALE(RemapALE, RemapALEService);
/*---------------------------------------------------------------------------*/
