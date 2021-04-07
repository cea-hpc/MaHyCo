// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "../MahycoModule.h"
/**
 *******************************************************************************
 * \file computeGradPhiFace1()
 * \brief phase 1 de projection : premiere etapes
 *  calcul des gradients aux faces verticales ou horizontales suivant le cas
 *  calcul des longueurs des faces verticales ou horizontales suivant le cas
 *  calcul des largeurs de cellule dans le sens verticales ou horizontales
 *suivant le cas \param \return m_grad_phi, LfLagrange (pas utilie), m_h_cell_lagrange
 *******************************************************************************
 */
void MahycoModule::
computeGradPhiFace(Integer idir, String name)  {
  debug() << " Entree dans computeGradPhiFace()";
  String ajout_interne = "_INTERNE";
  name = name + ajout_interne;
  m_h_cell_lagrange.fill(0.0);
  if (options()->ordreProjection > 1) {
    FaceGroup inner_dir_faces = mesh()->faceFamily()->findGroup(name);
    ENUMERATE_FACE(iface, inner_dir_faces) {
      Face face = *iface; 
      Cell cellb = face.backCell();
      Cell cellf = face.frontCell();
      for (Integer ivar = 0 ; ivar <  m_nb_vars_to_project ; ++ivar) {
        m_grad_phi_face[iface][ivar] = (m_phi_lagrange[cellf][ivar] - m_phi_lagrange[cellb][ivar]) 
                                    / m_deltax_lagrange[iface];
        //if (iface.localId() == 4195) info() << ivar << " m_grad_phi_face[iface][ivar] " << m_grad_phi_face[iface][ivar]
        //    << " varcellf " << m_phi_lagrange[cellf][ivar] << " varcellb " << m_phi_lagrange[cellb][ivar];
      }
      // somme des distances entre le milieu de la maille et le milieu de la face
      m_h_cell_lagrange[cellb] +=  (m_face_coord[iface] - m_cell_coord[cellb]).abs();
      m_h_cell_lagrange[cellf] +=  (m_face_coord[iface] - m_cell_coord[cellf]).abs();     
    }
  }
}
/**
 *******************************************************************************
 * \file computeGradPhi()
 * \brief phase de projection : seconde etape
 *        calcul du gradient aux mailles limites
 *        calcul des flux pente-borne arriere ou avant aux mailles
 *                 à partir du gradient precedent
 * \param
 * \return gradPhi1, deltaPhiFaceAr, deltaPhiFaceAv
 *******************************************************************************
 */
void MahycoModule::computeGradPhiCell(Integer idir) {
  debug() << " Entree dans computeGradPhiCell()";
  Real3 dirproj = {0.5 * (1-idir) * (2-idir), 
                   1.0 * idir * (2 -idir), 
                   -0.5 * idir * (1 - idir)};  
  if (options()->ordreProjection > 1) {
    ENUMERATE_CELL(icell,allCells()) {
      Cell cell = * icell;
      Cell backcell = cell; // pour les mailles de bord
      Cell frontcell = cell; // pour les mailles de bord
      Face backFace;
      Face frontFace;
      Integer cas_possible(0);
      ENUMERATE_FACE(iface, cell.faces()){
        const Face& face = *iface;
//         info() << " face " << face.localId() << " est en direction : " << idir << " : "  
//         << (Integer) m_is_dir_face[face][idir];
        if ( m_is_dir_face[face][idir] == true) {
          if (face.backCell() == cell) {
          frontFace = face;
          if (face.frontCell().localId() != -1) frontcell = face.frontCell();
          cas_possible++;
          } 
          else if ( face.frontCell() == cell) {
          backFace = face;
          if (face.backCell().localId() != -1) backcell = face.backCell();
          cas_possible++;
          }
        }
        if (cas_possible > 2) info() << " mauvais algo maille " << cell;
      }
      bool voisinage_pure = (options()->projectionLimiteurMixte == true &&
      m_est_mixte[cell] == 0 && m_est_mixte[frontcell] == 0 && m_est_mixte[backcell] == 0 &&
      m_est_pure[cell] == m_est_pure[frontcell] && m_est_pure[cell] == m_est_pure[backcell] );
      int limiter = options()->projectionLimiteurId;
      if ((options()->projectionPenteBorne == 1) && voisinage_pure)
      limiter = options()->projectionLimiteurPureId;
       
      // calcul de m_grad_phi[cell] 
      computeAndLimitGradPhi(limiter, frontFace, backFace, cell, frontcell, backcell);/*
      if (cell.localId() == 1251) {
        for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {    
          info() << " cell " << cell.localId() << " ivar " << ivar << " = " << m_grad_phi[cell][ivar];  
          info() << " frontface " << frontFace.localId() << " = " << m_grad_phi_face[frontFace][ivar];
          info() << " backFace " << backFace.localId() << " = " << m_grad_phi_face[backFace][ivar];
        }
      }*/
      //info() << " REmap " << idir << " cell " << cell.localId();
      //for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) 
      //  info() << "phi " << ivar <<  "  " << m_phi_lagrange[cell][ivar];
    
      if (options()->projectionPenteBorne == 1) {
        // if (cstmesh->cylindrical_mesh) exy = varlp->faceNormal(flFaces);
        Real Flux_sortant_ar = math::dot(m_outer_face_normal[cell][backFace.localId()], dirproj) * m_face_normal_velocity[backFace];
        // if (cstmesh->cylindrical_mesh) exy = varlp->faceNormal(frFaces);
//
        Real Flux_sortant_av = math::dot(m_outer_face_normal[cell][frontFace.localId()], dirproj) * m_face_normal_velocity[frontFace];
        
        Real flux_dual = 0.5 * (m_face_normal_velocity[backFace] + m_face_normal_velocity[frontFace]);
        
        Integer calcul_flux_dual(0);
        if (options()->calcul_flux_masse == 2) calcul_flux_dual = 1;
        
        RealUniqueArray delta_phi_face(m_nb_vars_to_project);
        RealUniqueArray dual_phi_flux(m_nb_vars_to_project);
        
        // calcul de m_delta_phi_face_ar et m_dual_phi_flux
        if (voisinage_pure)
            computeFluxPPPure(cell, frontcell, backcell, Flux_sortant_ar, 
                            m_deltat_n, 0, options()->threshold, 
                            options()->projectionPenteBorneComplet, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux);
        else
            computeFluxPP(cell, frontcell, backcell, Flux_sortant_ar, 
                            m_deltat_n, 0, options()->threshold, 
                            options()->projectionPenteBorneComplet, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux);
            
        for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {    
            m_delta_phi_face_ar[cell][ivar] = delta_phi_face[ivar];
            // et pour avoir un flux dual 2D
            m_dual_phi_flux[cell][ivar] = dual_phi_flux[ivar] * 
            0.5 * (m_face_length_lagrange[backFace][idir] + m_face_length_lagrange[frontFace][idir]);
        }
        // calcul de m_delta_phi_face_av
        if (voisinage_pure)
            computeFluxPPPure(cell, frontcell, backcell, Flux_sortant_av, 
                            m_deltat_n, 0, options()->threshold, 
                            options()->projectionPenteBorneComplet, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux);
        else
            computeFluxPP(cell, frontcell, backcell, Flux_sortant_av, 
                            m_deltat_n, 0, options()->threshold, 
                            options()->projectionPenteBorneComplet, flux_dual,
                            calcul_flux_dual, delta_phi_face, dual_phi_flux);
            
        for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {    
            m_delta_phi_face_av[cell][ivar] = delta_phi_face[ivar];
            // dual_phi_flux ne sert pas ici : le flux dual a dèjà été calculé
        }
      }
    }     
  }
}
/**
*****************************************************************************
 * \file computeUpwindFaceQuantitiesForProjection1()
 * \brief  phase 1 de projection : troisieme etape
 *        calcul de phiFace1
 *     qui contient la valeur reconstruite à l'ordre 1, 2 ou 3 des variables
 * projetees qui contient les flux des variables projetees avec l'option
 * pente-borne \param 
 * \return phiFace1
*******************************************************************************
*/
void MahycoModule::computeUpwindFaceQuantitiesForProjection(Integer idir, String name) {
  debug() << " Entree dans computeUpwindFaceQuantitiesForProjection()";
  String ajout_interne = "_INTERNE";
  name = name + ajout_interne;
  FaceGroup inner_dir_faces = mesh()->faceFamily()->findGroup(name);
  ICartesianMesh* m_cartesian_mesh = ICartesianMesh::getReference(mesh());
  CellDirectionMng cdm(m_cartesian_mesh->cellDirection(idir));
  ENUMERATE_FACE(iface, inner_dir_faces) {
      Face face = *iface; 
      Cell cellb = face.backCell();
      Cell cellf = face.frontCell();
    
      // phiFace1 correspond
      // à la valeur de phi(x) à la face pour l'ordre 2 sans plateau pente
      // à la valeur du flux (integration de phi(x)) pour l'ordre 2 avec
      // Plateau-Pente à la valeur du flux (integration de phi(x)) pour
      // l'ordre 3
      if (options()->ordreProjection <= 2) {
        if (options()->projectionPenteBorne == 0) {
          if (m_face_normal_velocity[face] * m_deltax_lagrange[face] > 0.0) {
              // phi_cb + (dot((x_f - x_cb), face_normal) * grad_phi_cb)   
            for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++)                                          
              m_phi_face[face][ivar] = m_phi_lagrange[cellb][ivar]
                    + math::dot((m_face_coord[face] - m_cell_coord[cellb]), 
                    m_face_normal[face]) * m_grad_phi[cellb][ivar];                                
          } else { 
            for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++)   
              m_phi_face[face][ivar] = m_phi_lagrange[cellf][ivar] 
                    + math::dot((m_face_coord[face] - m_cell_coord[cellf]), 
                    m_face_normal[face]) * m_grad_phi[cellf][ivar];
          }
        } else {
          for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {   
            m_phi_face[face][ivar] =
                m_delta_phi_face_av[cellb][ivar] - m_delta_phi_face_ar[cellf][ivar];
          }
        }
     } else if (options()->ordreProjection == 3) {      
         if (idir == 1) {
           // bug dans arcane pour la direction Y ?
           // pas de coherence entre backcell et ccb.previous()
           // il faut inverser cellb et cellf dans la direction Y
           Cell cell_tmp;
           cell_tmp = cellb;
           cellb = cellf;
           cellf = cell_tmp;
         } 
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
        Real vdt  = m_face_normal_velocity[face] * m_deltat_n;
        for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {   
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
        int nbmat = m_nb_env;   
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
}
/**
 *******************************************************************************
 * \file computeUremap1()
 * \brief phase 1 de projection : etape finale
 *        calcul de la variable Uremap1 par ajout ou retrait des flux
          Mises à jour de l'indicateur mailles mixtes
          calcul de la valeur de Phi issu de Uremap1
 * \param
 * \return Uremap1, varlp->mixte, varlp->pure, varlp->Phi
 *******************************************************************************
 */
void MahycoModule::computeUremap(Integer idir)  {
    debug() << " Entree dans computeUremap()";
    Real3 dirproj = {0.5 * (1-idir) * (2-idir), 
                   1.0 * idir * (2 -idir), 
                   -0.5 * idir * (1 - idir)};  
    int nbmat = m_nb_env;
    Real flux;
    ENUMERATE_CELL(icell,allCells()) {
      Cell cell = * icell;
      RealUniqueArray flux_face(m_nb_vars_to_project);
      flux_face.fill(0.);
      ENUMERATE_FACE(iface, cell.faces()){
        const Face& face = *iface;
        Integer i = iface.localId();
        
        // recuperation de la surface de la face 
        //m_face_length_lagrange[face][idir] 
 
        for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {   
            flux = computeRemapFlux(
                                options()->ordreProjection,
                                options()->projectionPenteBorne,
                                m_face_normal_velocity[face], m_face_normal[face],
                                m_face_length_lagrange[face][idir] , m_phi_face[face][ivar],
                                m_outer_face_normal[cell][i], dirproj, m_deltat_n);
            flux_face[ivar] += flux;
            
            // stockage des flux de masses (ivar = nbmat + imat, imat =0,..,nbmat) 
            // aux faces pour la quantite de mouvement de
            // Vnr
            // limitation à trois mat car m_flux_masse_face est real3
            if (ivar >= nbmat && ivar < 2*nbmat) m_flux_masse_face[cell][i][ivar] = flux;
            
            // doit etre FluxFace2(cCells, fFacesOfCellC) : m_flux_face[cell][i][ivar]
            
//           if (cdl->FluxBC > 0) {
//             // flux exterieur eventuel
//             //
//             reduction8 = reduction8 + (computeBoundaryFluxes(1, cCells, exy));
//             //
//           }
        }
 
      }
 
      for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {   
        m_u_lagrange[cell][ivar] -= flux_face[ivar];
      }

      // diagnostics et controle
      for (int imat = 0; imat < nbmat; imat++) {
        if (m_u_lagrange[cell][nbmat + imat] < 0.) {
          if (abs(m_u_lagrange[cell][nbmat + imat]) > 1.e2 * options()->threshold)
            info() << " cell " << cell.localId()
                    << " proj 1 --masse tres faiblement negative   "
                    << " soit " << m_u_lagrange[cell][nbmat + imat]
                    << " et volume " << m_u_lagrange[cell][imat];
          m_u_lagrange[cell][nbmat + imat] = 0.;
        }
        if (m_u_lagrange[cell][2 * nbmat + imat] < 0.) {
          if (abs(m_u_lagrange[cell][nbmat + imat]) > 1.e2 * options()->threshold)
            info() << " cell " << cell.localId()
                    << " --energie tres faiblement negative "
                    << " cell " << m_u_lagrange[cell][2 * nbmat + imat];
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
        if (options()->convnerg == 1)
          m_phi_lagrange[cell][3 * nbmat + 2] =
            m_u_lagrange[cell][3 * nbmat + 2] / somme_masse;

    } else {
        // somme_volume doit etre égale à m_cell_volume[cell]
      for (Integer ivar = 0; ivar < m_nb_vars_to_project; ivar++) {
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

