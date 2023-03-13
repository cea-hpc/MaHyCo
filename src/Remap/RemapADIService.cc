// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapADIService.h"

Integer RemapADIService::getOrdreProjection() { return options()->ordreProjection;}
bool RemapADIService::hasProjectionPenteBorne() { return options()->projectionPenteBorne;}
bool RemapADIService::hasProjectionSimplePente() { return options()->projectionSimplePente;}
bool RemapADIService::hasConservationEnergieTotale() { return options()->conservationEnergieTotale;}
bool RemapADIService::isEuler() {return options()->getIsEulerScheme();}
/**
 **************************************-*****************************************/
void RemapADIService::appliRemap(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) {
    
    synchronizeUremap();  
    synchronizeDualUremap();
    
    Integer idir(-1);
    m_cartesian_mesh = ICartesianMesh::getReference(mesh());
    
    for( Integer i=0; i< mesh()->dimension(); ++i){
      
      idir = (i + m_sens_projection())%(mesh()->dimension());
      // cas 2D : epaisseur de une maillage dans la direciton de projection
      if (m_cartesian_mesh->cellDirection(idir).globalNbCell() == 1) continue;
      // cas 1D : on debranche la projeciton suivant Y
      // if (idir == 1) continue;
      
      // info() << " projection direction " << idir;
      // calcul des gradients des quantites à projeter aux faces 
      computeGradPhiFace(idir, nb_vars_to_project, nb_env);
      // info() << " after computeGradPhiFace " << idir;
      // calcul des gradients des quantites à projeter aux cellules
      // (avec limiteur ordinaire) 
      // et pour le pente borne, calcul des flux aux faces des cellules
      computeGradPhiCell(idir, nb_vars_to_project, nb_env);
      // info() << " after  computeGradPhiCell " << idir;
      // calcul de m_phi_face
      // qui contient la valeur reconstruite à l'ordre 1, 2 ou 3 des variables projetees 
      // et qui contient les flux des variables projetees avec l'option pente-borne
      computeUpwindFaceQuantitiesForProjection(idir, nb_vars_to_project, nb_env);
      // info() << " after  computeUpwindFaceQuantitiesForProjection " << idir;
      
      
      computeUremap(idir, nb_vars_to_project, nb_env, withDualProjection);
      // info() << " after  computeUremap " << idir;
      synchronizeUremap();
      
      if (withDualProjection) {
        computeDualUremap(idir, nb_env);
        synchronizeDualUremap();
      }
    }
    m_sens_projection = m_sens_projection()+1;
    m_sens_projection = m_sens_projection()%(mesh()->dimension());
    
    // recuperation des quantités aux cells et aux envcell
    remapVariables(dimension,  withDualProjection,  nb_vars_to_project,  nb_env);
    // info() << " after  remapVariables " ;
    
}
/**
 *******************************************************************************/
void RemapADIService::resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env) {
    
    
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
void RemapADIService::computeGradPhiFace(Integer idir, Integer nb_vars_to_project, Integer nb_env)  {
  
  debug() << " Entree dans computeGradPhiFace()";
  m_h_cell_lagrange.fill(0.0);
  
  FaceDirectionMng fdm(m_cartesian_mesh->faceDirection(idir));
  ENUMERATE_FACE(iface, fdm.allFaces()) {
    Face face = *iface; 
    m_is_dir_face[face][idir] = true;
  }
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
      m_h_cell_lagrange[cellb] +=  (m_face_coord[iface] - m_cell_coord[cellb]).normL2();
      m_h_cell_lagrange[cellf] +=  (m_face_coord[iface] - m_cell_coord[cellf]).normL2();
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
void RemapADIService::computeGradPhiCell(Integer idir, Integer nb_vars_to_project, Integer nb_env) {
    
  debug() << " Entree dans computeGradPhiCell()";
  Real deltat = m_global_deltat();
  Real3 dirproj = {0.5 * (1-idir) * (2-idir), 
                   1.0 * idir * (2 -idir), 
                   -0.5 * idir * (1 - idir)};  
  m_delta_phi_face_av.fill(0.0);
  m_delta_phi_face_ar.fill(0.0);
  m_grad_phi.fill(0.0);
  FaceDirectionMng fdm(m_cartesian_mesh->faceDirection(idir));
  if (options()->ordreProjection > 1) {
    ENUMERATE_CELL(icell,allCells()) {
      Cell cell = * icell;
      Cell backcell = cell; // pour les mailles de bord
      Cell frontcell = cell; // pour les mailles de bord
      Face backFace;
      Face frontFace;
      Integer indexbackface;
      Integer indexfrontface;
      Integer cas_possible(0);
      ENUMERATE_FACE(iface, cell.faces()){
        const Face& face = *iface;
        
        if ( m_is_dir_face[face][idir] == true) {
          DirFace dir_face = fdm[face];
          if (dir_face.previousCell() == cell) {
          frontFace = face;
          indexfrontface = iface.index(); 
          if (dir_face.nextCell().localId() != -1) frontcell = dir_face.nextCell();
          cas_possible++;
          } 
          else if ( dir_face.nextCell() == cell) {
          backFace = face;
          indexbackface = iface.index(); 
          if (dir_face.previousCell().localId() != -1) backcell = dir_face.previousCell();
          cas_possible++;
          }
        }
        if (cas_possible > 2) info() << " mauvais algo maille " << cell.localId();
        if (cas_possible == 2) continue;
      }
      bool voisinage_pure = (options()->projectionPenteBorneMixte == true &&
        m_est_mixte[cell] == 0 && m_est_mixte[frontcell] == 0 && m_est_mixte[backcell] == 0 &&
        m_est_pure[cell] == m_est_pure[frontcell] && m_est_pure[cell] == m_est_pure[backcell] );
      
      int limiter = options()->projectionLimiteurId;
      if ((options()->getProjectionLimiteurPureId() == 1) && voisinage_pure)
        limiter = options()->projectionLimiteurPureId;
      // calcul de m_grad_phi[cell] 
      computeAndLimitGradPhi(limiter, frontFace, backFace, cell, frontcell, backcell, nb_vars_to_project);

      if (options()->projectionPenteBorne == 1) {
        // if (cstmesh->cylindrical_mesh) exy = varlp->faceNormal(flFaces);
        Real Flux_sortant_ar = math::dot(m_outer_face_normal[cell][indexbackface], dirproj) * m_face_normal_velocity[backFace];
        // if (cstmesh->cylindrical_mesh) exy = varlp->faceNormal(frFaces);
        Real Flux_sortant_av = math::dot(m_outer_face_normal[cell][indexfrontface], dirproj) * m_face_normal_velocity[frontFace];
        
        Real flux_dual = 0.5 * (Flux_sortant_ar + Flux_sortant_av);
        
        Integer calcul_flux_dual(1); // a supprimer - calcul du flux dual fait en pente borne
        
        RealUniqueArray delta_phi_face(nb_vars_to_project);
        RealUniqueArray dual_phi_flux(nb_vars_to_project);
       
        // calcul de m_delta_phi_face_ar et m_dual_phi_flux
        if (voisinage_pure)
            computeFluxPPPure(cell, frontcell, backcell, Flux_sortant_ar, 
                            deltat, 0, options()->threshold, 
                            options()->projectionPenteBorneDebarFix, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux,
                            nb_env, nb_vars_to_project);
        else
            computeFluxPP(cell, frontcell, backcell, Flux_sortant_ar, 
                            deltat, 0, options()->threshold, 
                            options()->projectionPenteBorneDebarFix, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux,
                            nb_env, nb_vars_to_project);
            
        for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {    
            m_delta_phi_face_ar[cell][ivar] = delta_phi_face[ivar];
            // et pour avoir un flux dual 2D
            m_dual_phi_flux[cell][ivar] = dual_phi_flux[ivar] * 
            0.5 * (m_face_length_lagrange[backFace][idir] + m_face_length_lagrange[frontFace][idir]);
        }
        // calcul de m_delta_phi_face_av
        if (voisinage_pure)
            computeFluxPPPure(cell, frontcell, backcell, Flux_sortant_av, 
                            deltat, 1, options()->threshold, 
                            options()->projectionPenteBorneDebarFix, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux,
                            nb_env, nb_vars_to_project);
        else
            computeFluxPP(cell, frontcell, backcell, Flux_sortant_av, 
                            deltat, 1, options()->threshold, 
                            options()->projectionPenteBorneDebarFix, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux,
                            nb_env, nb_vars_to_project);
            
        for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {    
            m_delta_phi_face_av[cell][ivar] = delta_phi_face[ivar];
            // dual_phi_flux ne sert pas ici : le flux dual a dèjà été calculé
        }
      }
    }     
  }
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
void RemapADIService::computeUpwindFaceQuantitiesForProjection(Integer idir, Integer nb_vars_to_project, Integer nb_env) {
    
  debug() << " Entree dans computeUpwindFaceQuantitiesForProjection()";
  Real deltat = m_global_deltat();
  CellDirectionMng cdm(m_cartesian_mesh->cellDirection(idir));
  FaceDirectionMng fdm(m_cartesian_mesh->faceDirection(idir));
  m_phi_face.fill(0.0);
  Integer order2 = options()->ordreProjection - 1;
  ENUMERATE_FACE(iface, fdm.innerFaces()) {
      Face face = *iface; 
      DirFace dir_face = fdm[face];
      Cell cellb = dir_face.previousCell();
      Cell cellf = dir_face.nextCell();
      // phiFace1 correspond
      // à la valeur de phi(x) à la face pour l'ordre 2 sans plateau pente
      // à la valeur du flux (integration de phi(x)) pour l'ordre 2 avec
      // Plateau-Pente à la valeur du flux (integration de phi(x)) pour
      // l'ordre 3
      if (options()->ordreProjection <= 2) {
        if (options()->projectionPenteBorne == 0) {
          if (m_face_normal_velocity[face] * m_deltax_lagrange[face] > 0.0) {
              // phi_cb + (dot((x_f - x_cb), face_normal) * grad_phi_cb)   
            for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++)                                          
              m_phi_face[face][ivar] = m_phi_lagrange[cellb][ivar]
                    + order2 * math::dot((m_face_coord[face] - m_cell_coord[cellb]), 
                    m_face_normal[face]) * m_grad_phi[cellb][ivar];                                
          } else { 
            for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++)   
              m_phi_face[face][ivar] = m_phi_lagrange[cellf][ivar] 
                    + order2 * math::dot((m_face_coord[face] - m_cell_coord[cellf]), 
                    m_face_normal[face]) * m_grad_phi[cellf][ivar];
          }
        } else {         
          for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {   
            m_phi_face[face][ivar] = (m_delta_phi_face_av[cellb][ivar] - m_delta_phi_face_ar[cellf][ivar]);
          }
        }
     } else if (options()->ordreProjection == 3) {    
        Cell cellbb = cellb;
        Cell cellbbb = cellb;
        Cell cellff = cellf;
        Cell cellfff = cellf;
        DirCell ccb(cdm.cell(cellb));
        if (ccb.previous().localId() != -1) {
          cellbb = ccb.previous();
          cellbbb = cellbb;
          DirCell ccbb(cdm.cell(cellbb));
          if (ccbb.previous().localId() != -1) {
            cellbbb = ccbb.previous();
          }
        }
        DirCell ccf(cdm.cell(cellf));
        if (ccf.next().localId() != -1) {
          cellff = ccf.next();
          cellfff = cellff;
          DirCell ccff(cdm.cell(cellff));
          if (ccff.next().localId() != -1) {
            cellfff = ccff.next();   
          }
        } 
        Real vdt  = m_face_normal_velocity[face] * deltat;
        for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {   
            m_phi_face[face][ivar] = ComputeFluxOrdre3(
                m_phi_lagrange[cellbbb][ivar],
                m_phi_lagrange[cellbb][ivar],
                m_phi_lagrange[cellb][ivar],
                m_phi_lagrange[cellf][ivar],
                m_phi_lagrange[cellff][ivar],
                m_phi_lagrange[cellfff][ivar],
                m_h_cell_lagrange[cellbbb],
                m_h_cell_lagrange[cellbb],
                m_h_cell_lagrange[cellb],
                m_h_cell_lagrange[cellf],
                m_h_cell_lagrange[cellff],
                m_h_cell_lagrange[cellfff],
                vdt);
        } 
        bool voisinage_pure = (m_est_mixte[cellb] == 0 && m_est_mixte[cellf] == 0 && 
        m_est_pure[cellb] == m_est_pure[cellf]);
        int nbmat = nb_env;   
        if (!voisinage_pure) {
          // comme dans le pente borne, on evite le pb de debar sur les maille a voisinage mixte
          // pinfo() << " voisinage mixte " << cellb.localId() << " et " << cellf.localId();
          for (int imat = 0; imat < nbmat; imat++) {
            if (vdt>0. && m_phi_lagrange[cellb][imat]!=0. && m_phi_lagrange[cellb][nbmat+imat]!=0.) {
             m_phi_face[face][nbmat+imat] = (m_phi_lagrange[cellb][nbmat+imat]/m_phi_lagrange[cellb][imat])
                *m_phi_face[face][imat];
             m_phi_face[face][2*nbmat+imat] = (m_phi_lagrange[cellb][2*nbmat+imat]/m_phi_lagrange[cellb][nbmat+imat])
                *m_phi_face[face][nbmat+imat];
            } else if (vdt<0. && m_phi_lagrange[cellf][imat]!=0. && m_phi_lagrange[cellf][nbmat+imat]!=0.) {
             m_phi_face[face][nbmat+imat] = (m_phi_lagrange[cellf][nbmat+imat]/m_phi_lagrange[cellf][imat])
                *m_phi_face[face][imat];
             m_phi_face[face][2*nbmat+imat] = (m_phi_lagrange[cellf][2*nbmat+imat]/m_phi_lagrange[cellf][nbmat+imat])
                *m_phi_face[face][nbmat+imat];
            }
          }
        }
     }
  }
  m_phi_face.synchronize();      
  debug() << " fin de computeUpwindFaceQuantitiesForProjection()";
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
void RemapADIService::computeUremap(Integer idir, Integer nb_vars_to_project, Integer nb_env, Integer withDualProjection) {
    
    debug() << " Entree dans computeUremap " << nb_vars_to_project << " " << nb_env;
    Real3 dirproj = {0.5 * (1-idir) * (2-idir), 
                   1.0 * idir * (2 -idir), 
                   -0.5 * idir * (1 - idir)};  
    int nbmat = nb_env;
    Real deltat = m_global_deltat();
    Real flux;
    m_dual_phi_flux.fill(0.);
    ENUMERATE_CELL(icell,allCells()) {
      Cell cell = * icell;
      RealUniqueArray flux_face(nb_vars_to_project);
      flux_face.fill(0.);
      ENUMERATE_FACE(iface, cell.faces()){
        const Face& face = *iface;
        Integer i = iface.index(); 
        if (std::fabs(math::dot(m_face_normal[face], dirproj)) >= 1.0E-10) {
          Real face_normal_velocity(m_face_normal_velocity[face]);
          Real face_length(m_face_length_lagrange[face][idir]);
          Real3 outer_face_normal(m_outer_face_normal[cell][i]);
          Real outer_face_normal_dir = math::dot(outer_face_normal, dirproj);
          if (options()->projectionPenteBorne == 0) {
            for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {  
                flux = outer_face_normal_dir * face_normal_velocity * face_length * deltat * m_phi_face[face][ivar];
                flux_face[ivar] += flux;
                m_dual_phi_flux[cell][ivar] += 0.5 * flux * outer_face_normal[idir] ;
            }
            
          } else { 
            for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {  
                flux = outer_face_normal_dir * face_length * m_phi_face[face][ivar];
                flux_face[ivar] += flux;
                m_dual_phi_flux[cell][ivar] += 0.5 * flux * outer_face_normal[idir];
            }
            
         }
        }
      }
 
      for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {   
        // flux dual comme demi somme des flux des deux faces contributives dans la direction idir
       // m_dual_phi_flux[cell][ivar] = 0.5 * flux_face[ivar]; 
        m_u_lagrange[cell][ivar] -= flux_face[ivar];
      }
      
      if (options()->projectionPenteBorne == 1 && withDualProjection) {
        // dans le cas du pente borne i.e. projection d'ordre 2 en temps
        // on peut obtenir les volumes, masses ou energies projetes tres faible 
        // on les tronque
        // info() << "rentrée ici que pour les cas avec lagrange : pas les advections dans le vide" ;
        for (int imat = 0; imat < nbmat; imat++) {
            if (m_u_lagrange[cell][imat] <  1. * options()->threshold) {
                m_u_lagrange[cell][imat] = 0.;
                m_u_lagrange[cell][nbmat + imat] = 0.;        
                m_u_lagrange[cell][2 * nbmat + imat] = 0.;
            } 
            if (m_u_lagrange[cell][nbmat + imat] <  1. * options()->threshold) {
                m_u_lagrange[cell][imat] = 0.;
                m_u_lagrange[cell][nbmat + imat] = 0.;        
                m_u_lagrange[cell][2 * nbmat + imat] = 0.;
            }
            if (m_u_lagrange[cell][2*nbmat + imat] <  1. * options()->threshold) {
                m_u_lagrange[cell][imat] = 0.;
                m_u_lagrange[cell][nbmat + imat] = 0.;        
                m_u_lagrange[cell][2 * nbmat + imat] = 0.;
            }
        }
      } 
      // diagnostics et controle
      for (int imat = 0; imat < nbmat; imat++) {
        if (m_u_lagrange[cell][nbmat + imat] < 0.) {
          if (abs(m_u_lagrange[cell][nbmat + imat]) > 1.e2 * options()->threshold)
            info() << " cell " << cell.localId()
                    << " proj 1 --masse tres faiblement negative   "
                    << " soit " << m_u_lagrange[cell][nbmat + imat]
                    << " et volume " << m_u_lagrange[cell][imat];
                    
          // m_u_lagrange[cell][imat] = 0.;
          m_u_lagrange[cell][nbmat + imat] = 0.;
        }
        if (m_u_lagrange[cell][2 * nbmat + imat] < 0.) {
          if (abs(m_u_lagrange[cell][nbmat + imat]) > 1.e2 * options()->threshold)
            info() << " cell " << cell.localId()
                    << " --energie tres faiblement negative "
                    << " cell " << m_u_lagrange[cell][2 * nbmat + imat];
          // m_u_lagrange[cell][imat] = 0.;
          // m_u_lagrange[cell][nbmat + imat] = 0.;
          m_u_lagrange[cell][2 * nbmat + imat] = 0.;
        }
      }
      // Calcul du volume de la maille apres 
      double somme_volume = 0.;
      for (int imat = 0; imat < nbmat; imat++) {
        somme_volume += m_u_lagrange[cell][imat];
      }
          
 
      if (options()->projectionPenteBorne == 1) {
        // option ou on ne regarde pas la variation de rho, V et e
        // phi = (f1, f2, rho1*f1, rho2*f2, Vx, Vy, e1, e2
        // ce qui permet d'ecrire le flux telque
        // Flux = (dv1 = f1dv, dv2=f2*dv, dm1=rho1*df1, dm2=rho2*df2, d(mVx) =
        // Vx*(dm1+dm2), d(mVy) = Vy*(dm1+dm2), d(m1e1) = e1*dm1,  d(m2e2) =
        // e2*dm2 dans computeFluxPP
        // Phi volume
        double somme_masse = 0.;
        for (int imat = 0; imat < nbmat; imat++) {
          m_phi_lagrange[cell][imat] = m_u_lagrange[cell][imat] / somme_volume;
        // Phi masse
        if (m_u_lagrange[cell][imat] != 0.)
          m_phi_lagrange[cell][nbmat + imat] =
                m_u_lagrange[cell][nbmat + imat] / m_u_lagrange[cell][imat];
        else
          m_phi_lagrange[cell][nbmat + imat] = 0.;
        somme_masse += m_u_lagrange[cell][nbmat + imat];
        }
        if (somme_masse!=0) {
          // Phi Vitesse
          m_phi_lagrange[cell][3 * nbmat] =
                m_u_lagrange[cell][3 * nbmat] / somme_masse;
          m_phi_lagrange[cell][3 * nbmat + 1] =
                m_u_lagrange[cell][3 * nbmat + 1] / somme_masse;
        }
        // Phi energie
        for (int imat = 0; imat < nbmat; imat++) {
        if (m_u_lagrange[cell][nbmat + imat] != 0.)
            m_phi_lagrange[cell][2 * nbmat + imat] =
                m_u_lagrange[cell][2 * nbmat + imat] /
                m_u_lagrange[cell][nbmat + imat];
        else
            m_phi_lagrange[cell][2 * nbmat + imat] = 0.;
        }
        // Phi energie cinétique
        if (options()->conservationEnergieTotale == 1)
          m_phi_lagrange[cell][3 * nbmat + 2] =
            m_u_lagrange[cell][3 * nbmat + 2] / somme_masse;

    } else {
        // somme_volume doit etre égale à m_cell_volume[cell]
      for (Integer ivar = 0; ivar < nb_vars_to_project; ivar++) {
        m_phi_lagrange[cell][ivar] = m_u_lagrange[cell][ivar] / somme_volume;
      }
    }
    // Mises à jour de l'indicateur mailles mixtes   
    Integer imatcell(0);
    Integer imatpure(-1);  
    for (int imat = 0; imat < nbmat; imat++)
      if (m_phi_lagrange[cell][imat] > 0.) {
        imatcell++;
        imatpure = imat;
      }  
    if (imatcell > 1) {
      m_est_mixte[cell] = 1;
      m_est_pure[cell] = -1;
    } else {
      m_est_mixte[cell] = 0;
      m_est_pure[cell] = imatpure;
    }
  }
}
/**
 *******************************************************************************
 * \file synchronizeUremap()
 * \brief phase de synchronisation des variables de projection apres projection
 * \return m_phi_lagrange, m_u_lagrange synchonise sur les mailles fantomes
 *******************************************************************************
 */
void RemapADIService::synchronizeUremap()  {
    debug() << " Entree dans synchronizeUremap()";
    m_phi_lagrange.synchronize();
    m_u_lagrange.synchronize();
    m_est_mixte.synchronize();
    m_est_pure.synchronize();
    m_dual_phi_flux.synchronize();
}/*---------------------------------------------------------------------------*/
ARCANE_REGISTER_SERVICE_REMAPADI(RemapADI, RemapADIService);
/*---------------------------------------------------------------------------*/
