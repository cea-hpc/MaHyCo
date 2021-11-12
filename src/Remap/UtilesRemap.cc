// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapADIService.h"

/**
 *******************************************************************************
 * \file computeDualVerticalGradPhi
 * \brief Calcul du grandient des variables duales 
 *        par appel à computeAndLimitGradPhi
 *        construction des grandients top et bottom et des largeurs h des 
 *        mailles duales
 * \param  pNodes
 * \return valeur du gradient dual et limité
 *******************************************************************************
 */
void RemapADIService::computeDualGradPhi(Node inode, 
                                          Node frontfrontnode, 
                                          Node frontnode, 
                                          Node backnode, 
                                          Node backbacknode, Integer idir) {
    
  Real3 grad_front(0. , 0. , 0.);
  Real3 grad_back(0. , 0. , 0.);
  int limiter = options()->projectionLimiteurId;
  
  // gradient vitesse X selon la direction idir
  grad_front[0] = (m_phi_dual_lagrange[frontnode][0] - m_phi_dual_lagrange[inode][0]) /
                (m_node_coord[frontnode][idir] - m_node_coord[inode][idir]);

  grad_back[0] = (m_phi_dual_lagrange[inode][0] - m_phi_dual_lagrange[backnode][0]) /
                (m_node_coord[inode][idir] - m_node_coord[backnode][idir]);

  // gradient vitesse Y selon la direction idir
  grad_front[1] = (m_phi_dual_lagrange[frontnode][1] - m_phi_dual_lagrange[inode][1]) /
                (m_node_coord[frontnode][idir] - m_node_coord[inode][idir]);

  grad_back[1] = (m_phi_dual_lagrange[inode][1] - m_phi_dual_lagrange[backnode][1]) /
                (m_node_coord[inode][idir] - m_node_coord[backnode][idir]);

  // gradient vitesse Z selon la direction idir
  grad_front[2] = (m_phi_dual_lagrange[frontnode][2] - m_phi_dual_lagrange[inode][2]) /
                (m_node_coord[frontnode][idir] - m_node_coord[inode][idir]);

  grad_back[2] = (m_phi_dual_lagrange[inode][2] - m_phi_dual_lagrange[backnode][2]) /
                (m_node_coord[inode][idir] - m_node_coord[backnode][idir]);

  // largeurs des mailles duales
  Real hmoins, h0, hplus;
  h0 = 0.5 * (m_node_coord[frontnode][idir]- m_node_coord[backnode][idir]);
  if (backbacknode.localId() == -1) {
    hmoins = 0.;
    hplus = 0.5 * (m_node_coord[frontfrontnode][idir]- m_node_coord[inode][idir]);
  } else if (frontfrontnode.localId() == -1) {
    hplus = 0.;
    hmoins = 0.5 * (m_node_coord[inode][idir] - m_node_coord[backbacknode][idir]);
  } else {
    hmoins = 0.5 * (m_node_coord[inode][idir] - m_node_coord[backbacknode][idir]);
    hplus = 0.5 * (m_node_coord[frontfrontnode][idir]- m_node_coord[inode][idir]);
  }

  // la fonction appelee ici
  computeAndLimitGradPhiDual(limiter, 
                             inode, frontnode, backnode,
                             grad_front, grad_back, 
                             h0, hplus, hmoins);
  // std::cout << " Nodes " << pNode << " gradV " << res << std::endl;
  // std::cout << " Nodes " << BottomBottomNode << " " << BottomNode << " " <<
  //   pNode << " " <<  TopNode << " " << TopTopNode << std::endl;   if (projectionLimiterId < minmodG) {
}


/**
 *******************************************************************************
 * \file INTY
 * \brief calcul de l'integrale de -infini à X de la fonction lineaire
 *    1_[x0,x1] ((x1 − x)y0 + (x − x0)y1)/(x1-x0)
 *    valant :
 *     y0 + 0.5 δ(y1 − y0)(x1 − x0)δ ou δ = S_[0,1]((X − x0)/(x1 − x0))
 *     zt S_[0,1](x) = min{max{x, 0}, 1}
 * \param  X, x0, y0, x1, y1
 * \return  valeur de l'integral
 *******************************************************************************
 */
 Real RemapADIService::INTY(Real X, Real x0, Real y0, Real x1, Real y1) {
  Real flux = 0.;
  // std::cout << " x0 " << x0 << std::endl;
  // std::cout << " x1 " << x1 << std::endl;
  if (abs(x1 - x0) > options()->threshold) {
    Real eta =
        math::min(math::max(0., (X - x0) / (x1 - x0)), 1.);
    // std::cout << " eta " << eta << std::endl;
    flux = (y0 + 0.5 * eta * (y1 - y0)) * (x1 - x0) * eta;
  }
  return flux;
}
/**
 *******************************************************************************
 * \file fluxLimiter
 * \brief retourne le valeur du limiteur demandé par la méthode classique
 *        qui suppose que les 3 mailles (utilise pour le calcul des gradient)
 *        ont la meme largeur
 * \param  r = gradmoins / gradplus ou gradplus / gradmoins
 * \return valeur du limiteur demandé
 *******************************************************************************
 */
 Real RemapADIService::fluxLimiter(int projectionLimiterId, Real r) {
  if (projectionLimiterId == minmod) {
    return std::max(0.0, std::min(1.0, r));
  } else if (projectionLimiterId == superBee) {
    return std::max(0.0, std::max(std::min(2.0 * r, 1.0), std::min(r, 2.0)));
  } else if (projectionLimiterId == vanLeer) {
    if (r <= 0.0)
      return 0.0;
    else
      return 2.0 * r / (1.0 + r);
  } else
    return 0.0;  // ordre 1
}
/**
 *******************************************************************************
 * \file fluxLimiterG
 * \brief retourne le valeur du limiteur demandé par la méthode exacte
 *        qui ne suppose pas que les 3 mailles (utilise pour le calcul des
 *gradient) ont la meme largeur (h0 != hplus != hmoins) \param  h0, hplus,
 *hmoins, gradplus, gradmoins \return valeur du limiteur demandé
 *******************************************************************************
 */
Real RemapADIService::fluxLimiterG(int projectionLimiterId, Real gradplus,
                           Real gradmoins, Real y0, Real yplus,
                           Real ymoins, Real h0, Real hplus,
                           Real hmoins) {
  Real grady, gradM, gradMplus, gradMmoins;
  // limitation rupture de pente (formule 16 si on utilise pas le plateau pente)
  if (gradplus * gradmoins < 0.0) return 0.;

  if (projectionLimiterId == minmodG)  // formule 9c
  {
    if ((yplus - ymoins) > 0.)
      grady = std::min(fabs(gradplus), fabs(gradmoins));
    else
      grady = -std::min(fabs(gradplus), fabs(gradmoins));
  } else if (projectionLimiterId == superBeeG)  // formule 9g
  {
    if ((yplus - ymoins) > 0.)
      grady = std::max(fabs(gradplus), fabs(gradmoins));
    else
      grady = -std::max(fabs(gradplus), fabs(gradmoins));
  } else if (projectionLimiterId == vanLeerG)  // formule 9e
  {
    Real lambdaplus = (h0 / 2. + hplus) / (h0 + hplus + hmoins);
    Real lambdamoins = (h0 / 2. + hmoins) / (h0 + hplus + hmoins);
    if ((lambdaplus * gradplus + lambdamoins * gradmoins) != 0.) {
      grady = gradplus * gradmoins /
              (lambdaplus * gradplus + lambdamoins * gradmoins);
    } else
      grady = 0.;
  } else if (projectionLimiterId == ultrabeeG) {
    grady = (yplus - ymoins) / h0;
  } else if (projectionLimiterId == arithmeticG) {
    Real lambdaplus = (h0 / 2. + hplus) / (h0 + hplus + hmoins);
    Real lambdamoins = (h0 / 2. + hmoins) / (h0 + hplus + hmoins);
    grady = lambdamoins * gradplus + lambdaplus * gradmoins;
  }
  // limitation simple-pente (formule 10)
  gradMplus = gradplus * (h0 + hplus) / h0;
  gradMmoins = gradmoins * (h0 + hmoins) / h0;
  gradM = std::min(fabs(gradMplus), fabs(gradMmoins));
  if ((yplus - ymoins) > 0.)
    grady = std::min(fabs(gradM), fabs(grady));
  else
    grady = -std::min(fabs(gradM), fabs(grady));

  return grady;
}
/**
 *******************************************************************************
 * \file computeY0
 * \brief calcul des Seuils de monotonie des reconstructions lineraires
 *  simple-pente pour l'option pente-borne
 *
 * \param  h0, hplus, hmoins, y0, yplus, ymoins
 * \return y0plus, y0moins
 *******************************************************************************
 */
Real RemapADIService::computeY0(int projectionLimiterId, Real y0, Real yplus,
                        Real ymoins, Real h0, Real hplus, Real hmoins,
                        int type) {
  // retourne {{y0plus, y0moins}}
  Real y0plus = 0., y0moins = 0.;
  if (projectionLimiterId == minmodG ||
      projectionLimiterId == minmod)  // minmod
  {
    y0plus = yplus;
    y0moins = ymoins;
  } else if (projectionLimiterId == superBeeG ||
             projectionLimiterId == superBee)  // superbee
  {
    y0plus = ((h0 + hmoins) * yplus + h0 * ymoins) / (2 * h0 + hmoins);
    y0moins = ((h0 + hplus) * ymoins + h0 * yplus) / (2 * h0 + hplus);
  } else if (projectionLimiterId == vanLeerG ||
             projectionLimiterId == vanLeer)  // vanleer
  {
    Real a = math::min(yplus, ymoins);
    Real b = math::max(yplus, ymoins);
    Real xplus = (h0 * h0 + 3 * h0 * hmoins + 2 * hmoins * hmoins) * yplus;
    Real xmoins = (h0 * h0 + 3 * h0 * hplus + 2 * hplus * hplus) * ymoins;
    xplus +=
        (h0 * h0 - h0 * hplus - 2 * hplus * hplus + 2 * h0 * hmoins) * ymoins;
    xmoins +=
        (h0 * h0 - h0 * hmoins - 2 * hmoins * hmoins + 2 * h0 * hplus) * yplus;
    xplus /= (2 * h0 * h0 + 5 * h0 * hmoins + 2 * hmoins * hmoins - h0 * hplus -
              2 * hplus * hplus);
    xmoins /= (2 * h0 * h0 + 5 * h0 * hplus + 2 * hplus * hplus - h0 * hmoins -
               2 * hmoins * hmoins);

    y0plus = math::min(math::max(xplus, a), b);
    y0moins = math::min(math::max(xmoins, a), b);
  } else if (projectionLimiterId == ultrabeeG) {
    y0plus = (yplus + ymoins) / 2.;
    y0moins = (yplus + ymoins) / 2.;
  } else if (projectionLimiterId == arithmeticG) {
    y0plus = ((h0 + hmoins + hplus) * yplus + h0 * ymoins) /
             (2 * h0 + hmoins + hplus);
    y0moins = ((h0 + hmoins + hplus) * ymoins + h0 * yplus) /
              (2 * h0 + hmoins + hplus);
  } else if (projectionLimiterId == 3000) {
    y0plus = yplus;
    y0moins = ymoins;
  }
  if (type == 0)
    return y0plus;
  else if (type == 1)
    return y0moins;
  else
    return 0.0;  // lancer forcement avec type 0 ou 1 mais warning compile
}
/**
 *******************************************************************************
 * \file computexgxd
 * \brief calcul des absisses des points d'appui de la reconstruction en 3
 *morceaux
 *
 * \param  h0, hplus, hmoins, y0, yplus, ymoins
 * \return xg, xd
 *******************************************************************************
 */
Real RemapADIService::computexgxd(Real y0, Real yplus, 
                                          Real ymoins, Real h0,
                          Real y0plus, Real y0moins, int type) {
  // retourne {{xg, xd}}
  Real xd = 0., xg = 0.;
  Real xplus = 1.;
  if (abs(y0plus - yplus) > options()->threshold)
    xplus = (y0 - yplus) / (y0plus - yplus) - 1. / 2.;
  Real xmoins = 1.;
  if (abs(y0moins - ymoins) > options()->threshold)
    xmoins = (y0 - ymoins) / (y0moins - ymoins) - 1. / 2.;
  xd = +h0 * math::min(math::max(xplus, -1. / 2.), 1. / 2.);
  xg = -h0 * math::min(math::max(xmoins, -1. / 2.), 1. / 2.);
  if (type == 0)
    return xg;
  else if (type == 1)
    return xd;
  else
    return 0.0;  // lancer forcement avec type 0 ou 1 mais warning compile
}
/**
 *******************************************************************************
 * \file computeygyd
 * \brief calcul des ordonnees des points d'appui de la reconstruction en 3
 *morceaux
 *
 * \param  h0, hplus, hmoins, y0, yplus, ymoins
 * \return yg, yd
 *******************************************************************************
 */
Real RemapADIService::computeygyd(Real y0, Real yplus, Real ymoins, Real h0,
                          Real y0plus, Real y0moins, Real grady,
                          int type) {
  // retourne {{yg, yd}}
  Real yd, yg;
  Real xtd = y0 + h0 / 2 * grady;
  Real xtg = y0 - h0 / 2 * grady;
  Real ad = math::min(yplus, 2. * y0moins - ymoins);
  Real bd = math::max(yplus, 2. * y0moins - ymoins);
  Real ag = math::min(ymoins, 2. * y0plus - yplus);
  Real bg = math::max(ymoins, 2. * y0plus - yplus);
  yd = math::min(math::max(xtd, ad), bd);
  yg = math::min(math::max(xtg, ag), bg);
  if (type == 0)
    return yg;
  else if (type == 1)
    return yd;
  else
    return 0.0;  // lancer forcement avec type 0 ou 1 mais warning compile
}
/**
 *******************************************************************************
 * \file computeAndLimitGradPhi
 * \brief Calcul du grandient des variables  (phi) à partir
 *        des mailles et des faces amonts : backcell et backFace
 *        et avals : frontcell et frontFace 
 * \param  
 * \return valeur du gradient limité de chaque variables phi
 *******************************************************************************
 */
void RemapADIService::
computeAndLimitGradPhi(int projectionLimiterId, Face frontFace, Face backFace, 
                       Cell cell, Cell frontcell, Cell backcell, int nb_vars) {
    
   if (projectionLimiterId < minmodG) {
     // info() << " Passage gradient limite Classique ";
    for (Integer ivar = 0; ivar < nb_vars; ivar++) {
      m_grad_phi[cell][ivar] = 0.;
      if (m_grad_phi_face[backFace][ivar] != 0.) 
        m_grad_phi[cell][ivar] += 0.5 * (
                fluxLimiter(projectionLimiterId, 
                    m_grad_phi_face[frontFace][ivar] /
                    m_grad_phi_face[backFace][ivar]) 
                    * m_grad_phi_face[backFace][ivar]);
     if (m_grad_phi_face[frontFace][ivar] !=0.) 
        m_grad_phi[cell][ivar] += 0.5 * (
                fluxLimiter(projectionLimiterId, 
                    m_grad_phi_face[backFace][ivar] /
                    m_grad_phi_face[frontFace][ivar]) 
                    * m_grad_phi_face[frontFace][ivar]);
    }
  } else {
    // info() << " Passage gradient limite Genéralisé ";
    for (Integer ivar = 0; ivar < nb_vars; ivar++) {
      m_grad_phi[cell][ivar] =
          fluxLimiterG(projectionLimiterId, m_grad_phi_face[frontFace][ivar], 
                       m_grad_phi_face[backFace][ivar], m_phi_lagrange[cell][ivar],
                       m_phi_lagrange[frontcell][ivar], m_phi_lagrange[backcell][ivar],
                       m_h_cell_lagrange[cell],
                       m_h_cell_lagrange[frontcell], m_h_cell_lagrange[backcell]);
    }
  }
}
/**
 ******************************************************************************** 
 * \file computeAndLimitGradPhiDual
 * \brief Calcul du grandient des variables (m_phi_dual_lagrange) à partir
 *        des noeuds amont : backcnode et aval : frontnode 
 * \param  grad_back et grad_front - gradients en amon et aval du noeud
 *          h0,  hplus,  hmoins : largeur des cellules duales
 * \return valeur du gradient limité  m_dual_grad_phi
 
 *******************************************************************************
 */
void RemapADIService::
computeAndLimitGradPhiDual(int projectionLimiterId, Node inode, 
                           Node frontnode, Node backnode, 
                           Real3 grad_front, Real3 grad_back, Real h0, Real hplus, Real hmoins) {
    
   if (projectionLimiterId < minmodG) {
     // info() << " Passage gradient limite Classique ";
    for (Integer ivar = 0; ivar < 3; ivar++) {
      m_dual_grad_phi[inode][ivar] = 0.;
      if (grad_back[ivar] != 0. ) 
        m_dual_grad_phi[inode][ivar] += 0.5 *
                fluxLimiter(projectionLimiterId, 
                    grad_front[ivar] / grad_back[ivar]) 
                    * grad_back[ivar];
      if (grad_front[ivar] !=0.)
         m_dual_grad_phi[inode][ivar] += 0.5 *
                fluxLimiter(projectionLimiterId, 
                    grad_back[ivar] / grad_front[ivar]) 
                    * grad_front[ivar];
    }
  } else {
    // info() << " Passage gradient limite Genéralisé ";
    for (Integer ivar = 0; ivar < 3; ivar++) {
      m_dual_grad_phi[inode][ivar] =
          fluxLimiterG(projectionLimiterId, grad_front[ivar], 
                       grad_back[ivar], 
                       m_phi_dual_lagrange[inode][ivar],
                       m_phi_dual_lagrange[frontnode][ivar], m_phi_dual_lagrange[backnode][ivar],
                       h0, hplus, hmoins);
    }
  }
}
/**
 *******************************************************************************
 * \file computeFluxPP
 * \brief Calcul des flux pente-borne à partir des valeurs des variables 
 *  sur 3 mailles moins, 0, plus
 *  On integre sur -h0/2 et -h0/2.+partie_positive_v pour le flux à gauche (type 0)
 *  On integre sur  h0/2.-partie_positive_v et h0/2 pour le flux à droite (type 1)
 *
 *  pour les mailles mixtes ou celles à voisinages mixtes 
 *     on deduit des flux des volumes partiels, les flux des masses partielles
 *     dm = rho * dv
 *     on deduit des flux des masses partielles, les flux d'energie et de vitesse
 *     d(me) = e * (dm)
 *     rho est soit la valeur dans la maille de la densite 
 *     soit et une valeur reconstruite a par des flux (option complet)
 *     e celle dans la maille de l'energie interne
 *     soit et une valeur reconstruite a par des flux (option complet)
 *
 *  On calcule de la meme facon des flux "duales" en integrant
 *  entre 0 et partie_positive_v pour le flux duale 
 * \param 
 * \return flux des variables phi
 *******************************************************************************
 */
void RemapADIService::computeFluxPP(Cell cell, Cell frontcell, Cell backcell, 
                                     Real face_normal_velocity, 
                                     Real deltat_n, Integer type, Real flux_threshold, 
                                     Integer projectionPenteBorneDebarFix, 
                                     Real dual_normal_velocity,
                                     Integer calcul_flux_dual,
                                     RealArrayView Flux, RealArrayView Flux_dual,
                                     int nbmat, int nb_vars
                                    ) {
    

  Flux.fill(0.0);
  Flux_dual.fill(0.0);

  Real y0plus, y0moins, xd, xg, yd, yg;
  Real flux1, flux2, flux3, flux1m, flux2m, flux3m;
  Real partie_positive_v = 0.5 * (face_normal_velocity + abs(face_normal_velocity)) * deltat_n;
  Real partie_positive_dual_v = 0.5 * (dual_normal_velocity + abs(dual_normal_velocity)) * deltat_n;
  Integer cas_PP = 0;
    for (Integer ivar = 0; ivar < nb_vars; ivar++) {
      Real h0 = m_h_cell_lagrange[cell];
      Real hplus = m_h_cell_lagrange[frontcell];
      Real hmoins = m_h_cell_lagrange[backcell];
      // calcul des seuils y0plus, y0moins pour cCells
      y0plus = computeY0(options()->projectionLimiteurId, 
                         m_phi_lagrange[cell][ivar],        
                         m_phi_lagrange[frontcell][ivar], 
                         m_phi_lagrange[backcell][ivar], 
                         h0, hplus, hmoins, 0);
      
      y0moins = computeY0(options()->projectionLimiteurId, 
                         m_phi_lagrange[cell][ivar],        
                         m_phi_lagrange[frontcell][ivar], 
                         m_phi_lagrange[backcell][ivar], 
                         h0, hplus, hmoins, 1);
      // calcul des points d'intersections xd,xg
      xg = computexgxd(m_phi_lagrange[cell][ivar], 
                        m_phi_lagrange[frontcell][ivar], 
                        m_phi_lagrange[backcell][ivar], 
                        h0, y0plus, y0moins, 0);
      xd = computexgxd(m_phi_lagrange[cell][ivar], 
                        m_phi_lagrange[frontcell][ivar], 
                        m_phi_lagrange[backcell][ivar], 
                        h0, y0plus, y0moins, 1);
      // calcul des valeurs sur ces points d'intersections
      yg = computeygyd(m_phi_lagrange[cell][ivar], 
                        m_phi_lagrange[frontcell][ivar], 
                        m_phi_lagrange[backcell][ivar], 
                        h0, y0plus, y0moins, m_grad_phi[cell][ivar], 0);
      yd = computeygyd(m_phi_lagrange[cell][ivar], 
                        m_phi_lagrange[frontcell][ivar], 
                        m_phi_lagrange[backcell][ivar], 
                        h0, y0plus, y0moins, m_grad_phi[cell][ivar], 1);
        //  
     if (type == 0)  
         // flux arriere ou en dessous de cCells, integration entre
         // -h0/2. et -h0/2.+abs(face_normal_velocity)*deltat_n
      {
      // Flux1m : integrale -inf,  -h0/2.+partie_positive_v
      flux1m =
          INTY(-h0 / 2. + partie_positive_v, -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      // Flux1m : integrale -inf,  -h0/2.
      flux1 = INTY(-h0 / 2., -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      // Flux2m : integrale -inf,  -h0/2.+partie_positive_v
      flux2m = INTY(-h0 / 2. + partie_positive_v, xg, yg, xd, yd);
      // Flux2 : integrale -inf,  -h0/2.
      flux2 = INTY(-h0 / 2., xg, yg, xd, yd);
      // Flux3m : integrale -inf,  -h0/2.+partie_positive_v
      flux3m = INTY(-h0 / 2. + partie_positive_v, xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
      // Flux3 : integrale -inf,  -h0/2.
      flux3 = INTY(-h0 / 2., xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
      // integrale positive
      Flux[ivar] = std::max(
          ((flux1m - flux1) + (flux2m - flux2) + (flux3m - flux3)), 0.);
      // formule 16
      if (((m_phi_lagrange[frontcell][ivar] - m_phi_lagrange[cell][ivar]) * 
          (m_phi_lagrange[backcell][ivar] - m_phi_lagrange[cell][ivar])) >= 0.)
        Flux[ivar] = m_phi_lagrange[cell][ivar] * partie_positive_v;
        
      //
      // et calcul du flux dual si calcul_flux_dual=1
      if (calcul_flux_dual == 1) {
        // Flux1m : integrale -inf, partie_positive_dual_v
        flux1m = INTY(partie_positive_dual_v, -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
        // Flux1m : integrale -inf,  0..
        flux1 = INTY(0., -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
        // Flux2m : integrale -inf, partie_positive_dual_v
        flux2m = INTY(partie_positive_dual_v, xg, yg, xd, yd);
        // Flux2 : integrale -inf, 0.
        flux2 = INTY(0., xg, yg, xd, yd);
        // Flux3m : integrale -inf, partie_positive_dual_v
        flux3m = INTY(partie_positive_dual_v, xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
        // Flux3 : integrale -inf, 0.
        flux3 = INTY(0., xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
        // integrale positive
        Flux_dual[ivar] = std::max(
            ((flux1m - flux1) + (flux2m - flux2) + (flux3m - flux3)), 0.);
        // formule 16
        if (((m_phi_lagrange[frontcell][ivar] - m_phi_lagrange[cell][ivar]) * 
          (m_phi_lagrange[backcell][ivar] - m_phi_lagrange[cell][ivar])) >= 0.)
          Flux_dual[ivar] = m_phi_lagrange[cell][ivar] * partie_positive_dual_v;
        //
      }
    } else if (type == 1) {
      // flux devant ou au dessus de cCells, integration entre
      // h0/2.-abs(face_normal_velocity)*deltat_n et h0/2.
      // Flux1 : integrale -inf,  h0/2.-partie_positive_v
      flux1 = INTY(h0 / 2. - partie_positive_v, -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      // Flux1m : integrale -inf,  -h0/2.
      flux1m = INTY(h0 / 2., -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      //
      // Flux2 : integrale -inf,  h0/2.-partie_positive_v
      flux2 = INTY(h0 / 2. - partie_positive_v, xg, yg, xd, yd);
      // Flux2m : integrale -inf,  -h0/2.
      flux2m = INTY(h0 / 2., xg, yg, xd, yd);
      //
      // Flux3 : integrale -inf,  h0/2.-partie_positive_v
      flux3 = INTY(h0 / 2. - partie_positive_v, xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar] );
      // Flux3m : integrale -inf,  -h0/2.
      flux3m = INTY(h0 / 2., xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar] );
      //
      // integrale positive
      Flux[ivar] = std::max(
          ((flux1m - flux1) + (flux2m - flux2) + (flux3m - flux3)), 0.);
      // formule 16
        if (((m_phi_lagrange[frontcell][ivar] - m_phi_lagrange[cell][ivar]) * 
          (m_phi_lagrange[backcell][ivar] - m_phi_lagrange[cell][ivar])) >= 0.)
        Flux[ivar] = m_phi_lagrange[cell][ivar] * partie_positive_v;
      //
      // et calcul du flux dual si calcul_flux_dual=1
      if (calcul_flux_dual == 1) {
        // flux dual deja calculé lors du premier appel à la fonction
        Flux_dual[ivar] = 0.;
      }
    }
  }
  if (projectionPenteBorneDebarFix == 1) {
    // pinfo() << "projectionPenteBorneDebarFix" << projectionPenteBorneDebarFix;
    // les flux de masse se déduisent des flux de volume en utilisant une valeur
    // moyenne de Rho calculée par le flux de masse / flux de volume Celles des
    // energies massiques avec une valeur moyenne de e calculée par le flux
    // d'energie / flux de volume Celles de quantité de mouvement avec une
    // valeur moyenne de u calculée par le flux de vitesse / flux de volume
    Real somme_flux_masse = 0.;
    Real somme_flux_volume = 0.;
    for (size_t imat = 0; imat < nbmat; imat++) somme_flux_volume += Flux[imat];

    if (std::fabs(somme_flux_volume) > flux_threshold ) {
      for (size_t imat = 0; imat < nbmat; imat++) {
        Flux[nbmat + imat] =
            (Flux[nbmat + imat] / somme_flux_volume) * Flux[imat];
        Flux[2. * nbmat + imat] =
            (Flux[2. * nbmat + imat] / somme_flux_volume) * Flux[nbmat + imat];
        somme_flux_masse += Flux[nbmat + imat];
      }

      Flux[3 * nbmat] = (Flux[3 * nbmat] / somme_flux_volume) *
                        somme_flux_masse;  // flux de quantité de mouvement x
      Flux[3 * nbmat + 1] =
          (Flux[3 * nbmat + 1] / somme_flux_volume) *
          somme_flux_masse;  // flux de quantité de mouvement y
      Flux[3 * nbmat + 2] =
          (Flux[3 * nbmat + 2] / somme_flux_volume) * somme_flux_masse;
      Flux[3 * nbmat + 3] =
          m_phi_lagrange[cell][3 * nbmat + 3] * somme_flux_volume;  // flux pour la pseudo VNR
    } else {
      Flux.fill(0.);
    }
  } else if (projectionPenteBorneDebarFix == 2) {
    // pinfo() << "projectionPenteBorneDebarFix" << projectionPenteBorneDebarFix;
    // les flux de masse, de quantité de mouvement et d'energie massique se
    // deduisent des flux de volumes avec la valeur de rho, e et u à la maille
    Real somme_flux_masse = 0.;
    Real somme_flux_volume = 0.;
    for (size_t imat = 0; imat < nbmat; imat++) {
      Flux[nbmat + imat] =
          m_phi_lagrange[cell][nbmat + imat] * Flux[imat];  // flux de masse de imat
      Flux[2 * nbmat + imat] =
          m_phi_lagrange[cell][2 * nbmat + imat] *
          Flux[nbmat + imat];  // flux de masse energy de imat
      somme_flux_masse += Flux[nbmat + imat];
      somme_flux_volume += Flux[imat];
    }
    Flux[3 * nbmat] =
        m_phi_lagrange[cell][3 * nbmat] * somme_flux_masse;  // flux de quantité de mouvement x
    Flux[3 * nbmat + 1] = m_phi_lagrange[cell][3 * nbmat + 1] *
                          somme_flux_masse;  // flux de quantité de mouvement y
    Flux[3 * nbmat + 2] =
        m_phi_lagrange[cell][3 * nbmat + 2] * somme_flux_masse;  // flux d'energie cinetique
    Flux[3 * nbmat + 3] =
        m_phi_lagrange[cell][3 * nbmat + 3] * somme_flux_volume;  // flux pour la pseudo VNR
  }

  if (partie_positive_v == 0.) Flux.fill(0.);
  if (partie_positive_dual_v == 0.) Flux_dual.fill(0.);

  return;
}
/**
 *******************************************************************************
 * \file computeFluxPPPure
 * \brief Calcul des flux pente-borne à partir des valeurs des variables 
 *  sur 3 mailles moins, 0, plus
 *  utilises pour les mailles pures à voisinages pures, 
 *  On integre sur -h0/2 et -h0/2.+partie_positive_v pour le flux à gauche (type 0)
 *  On integre sur  h0/2.-partie_positive_v et h0/2 pour le flux à droite (type 1)
 *
 *     on deduit des flux des volumes partiels, les flux des masses partielles
 *     dm = rho * dv
 *     on deduit des flux des masses partielles, les flux d'energie et de vitesse
 *     d(me) = e * (dm)
 *     rho est soit la valeur dans la maille de la densite 
 *     soit et une valeur reconstruite a par des flux (option complet)
 *     e celle dans la maille de l'energie interne
 *     soit et une valeur reconstruite a par des flux (option complet)
 *
 *  On calcule de la meme facon des flux "duales" en integrant
 *  entre 0 et partie_positive_v pour le flux duale 
 * \param 
 * \return Flux, Flux_dual
 *******************************************************************************
 */
void RemapADIService::computeFluxPPPure(Cell cell, Cell frontcell, Cell backcell, 
                                     Real face_normal_velocity, 
                                     Real deltat_n, Integer type, Real flux_threshold, 
                                     Integer projectionPenteBorneDebarFix, 
                                     Real dual_normal_velocity,
                                     Integer calcul_flux_dual,
                                     RealArrayView Flux, RealArrayView Flux_dual,
                                     int nbmat, int nb_vars
                                    ) {
  Flux.fill(0.0);
  Flux_dual.fill(0.0);   
  Real y0plus, y0moins, xd, xg, yd, yg;
  Real flux1, flux2, flux3, flux1m, flux2m, flux3m;
  Real partie_positive_v = 0.5 * (face_normal_velocity + abs(face_normal_velocity)) * deltat_n;
  Real partie_positive_dual_v = 0.5 * (dual_normal_velocity + abs(dual_normal_velocity)) * deltat_n;
  int cas_PP = 0;
  // on ne fait que la projection des volumes et masses
  for (Integer ivar = 0; ivar < nb_vars; ivar++) {
    Real h0 = m_h_cell_lagrange[cell];
    Real hplus = m_h_cell_lagrange[frontcell];
    Real hmoins = m_h_cell_lagrange[backcell];
    // calcul des seuils y0plus, y0moins pour cCells
    y0plus = computeY0(options()->projectionLimiteurPureId,
                         m_phi_lagrange[cell][ivar],        
                         m_phi_lagrange[frontcell][ivar], 
                         m_phi_lagrange[backcell][ivar], 
                         h0, hplus, hmoins, 0);
      
    y0moins = computeY0(options()->projectionLimiteurPureId,
                         m_phi_lagrange[cell][ivar],        
                         m_phi_lagrange[frontcell][ivar], 
                         m_phi_lagrange[backcell][ivar], 
                         h0, hplus, hmoins, 1);
    // calcul des points d'intersections xd,xg
    xg = computexgxd(m_phi_lagrange[cell][ivar], 
                      m_phi_lagrange[frontcell][ivar], 
                      m_phi_lagrange[backcell][ivar], 
                      h0, y0plus, y0moins, 0);
    xd = computexgxd(m_phi_lagrange[cell][ivar], 
                      m_phi_lagrange[frontcell][ivar], 
                      m_phi_lagrange[backcell][ivar], 
                      h0, y0plus, y0moins, 1);
    // calcul des valeurs sur ces points d'intersections
    yg = computeygyd(m_phi_lagrange[cell][ivar], 
                      m_phi_lagrange[frontcell][ivar], 
                      m_phi_lagrange[backcell][ivar], 
                      h0, y0plus, y0moins, m_grad_phi[cell][ivar], 0);
    yd = computeygyd(m_phi_lagrange[cell][ivar], 
                      m_phi_lagrange[frontcell][ivar], 
                      m_phi_lagrange[backcell][ivar], 
                      h0, y0plus, y0moins, m_grad_phi[cell][ivar], 1);
     //  

    if (type == 0)  // flux arriere ou en dessous de cCells, integration entre
                    // -h0/2. et -h0/2.+abs(face_normal_velocity)*deltat_n
    {
      // Flux1m : integrale -inf,  -h0/2.+partie_positive_v
      flux1m =
          INTY(-h0 / 2. + partie_positive_v, -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      // Flux1m : integrale -inf,  -h0/2.
      flux1 = INTY(-h0 / 2., -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      //
      // Flux2m : integrale -inf,  -h0/2.+partie_positive_v
      flux2m = INTY(-h0 / 2. + partie_positive_v, xg, yg, xd, yd);
      // Flux2 : integrale -inf,  -h0/2.
      flux2 = INTY(-h0 / 2., xg, yg, xd, yd);
      //
      // Flux3m : integrale -inf,  -h0/2.+partie_positive_v
      flux3m = INTY(-h0 / 2. + partie_positive_v, xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
      // Flux3 : integrale -inf,  -h0/2.
      flux3 = INTY(-h0 / 2., xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
      //
      // integrale positive
      Flux[ivar] = std::max(
          ((flux1m - flux1) + (flux2m - flux2) + (flux3m - flux3)), 0.);
      // formule 16
      if (((m_phi_lagrange[frontcell][ivar] - m_phi_lagrange[cell][ivar]) * 
          (m_phi_lagrange[backcell][ivar] - m_phi_lagrange[cell][ivar])) >= 0.) 
        Flux[ivar] = m_phi_lagrange[cell][ivar]  * partie_positive_v;
      //
      // et calcul du flux dual si calcul_flux_dual=1
      if (calcul_flux_dual == 1) {
        // Flux1m : integrale -inf, partie_positive_dual_v
        flux1m = INTY(partie_positive_dual_v, -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
        // Flux1m : integrale -inf,  0..
        flux1 = INTY(0., -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
        // Flux2m : integrale -inf, partie_positive_dual_v
        flux2m = INTY(partie_positive_dual_v, xg, yg, xd, yd);
        // Flux2 : integrale -inf, 0.
        flux2 = INTY(0., xg, yg, xd, yd);
        // Flux3m : integrale -inf, partie_positive_dual_v
        flux3m = INTY(partie_positive_dual_v, xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
        // Flux3 : integrale -inf, 0.
        flux3 = INTY(0., xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
        // integrale positive
        Flux_dual[ivar] = std::max(
            ((flux1m - flux1) + (flux2m - flux2) + (flux3m - flux3)), 0.);
        // formule 16
        if (((m_phi_lagrange[frontcell][ivar] - m_phi_lagrange[cell][ivar]) * 
          (m_phi_lagrange[backcell][ivar] - m_phi_lagrange[cell][ivar])) >= 0.)
          Flux_dual[ivar] = m_phi_lagrange[cell][ivar] * partie_positive_dual_v;
        //
      }
    } else if (type == 1) {
      // flux devant ou au dessus de cCells, integration entre
      // h0/2.-abs(face_normal_velocity)*deltat_n et h0/2.
      // Flux1 : integrale -inf,  h0/2.-partie_positive_v
      flux1 = INTY(h0 / 2. - partie_positive_v, -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      // Flux1m : integrale -inf,  -h0/2.
      flux1m = INTY(h0 / 2., -h0 / 2., m_phi_lagrange[backcell][ivar], xg, yg);
      // Flux2 : integrale -inf,  h0/2.-partie_positive_v
      flux2 = INTY(h0 / 2. - partie_positive_v, xg, yg, xd, yd);
      // Flux2m : integrale -inf,  -h0/2.
      flux2m = INTY(h0 / 2., xg, yg, xd, yd);
      //
      // Flux3 : integrale -inf,  h0/2.-partie_positive_v
      flux3 = INTY(h0 / 2. - partie_positive_v, xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
      // Flux3m : integrale -inf,  -h0/2.
      flux3m = INTY(h0 / 2., xd, yd, h0 / 2., m_phi_lagrange[frontcell][ivar]);
      // integrale positive
      Flux[ivar] = std::max(
          ((flux1m - flux1) + (flux2m - flux2) + (flux3m - flux3)), 0.);
      // formule 16
      if (((m_phi_lagrange[frontcell][ivar] - m_phi_lagrange[cell][ivar]) * 
          (m_phi_lagrange[backcell][ivar] - m_phi_lagrange[cell][ivar])) >= 0.)
        Flux[ivar] = m_phi_lagrange[cell][ivar] * partie_positive_v;
      // et calcul du flux dual si calcul_flux_dual=1
      if (calcul_flux_dual == 1) {
        // flux dual deja calculé lors du premier appel à la fonction
        Flux_dual[ivar] = 0.;
      }
    }
  }
  if (projectionPenteBorneDebarFix == 1) {
    // les flux d' energie se déduisent des flux de masse en utilisant 
    // une valeur moyenne de e calculée par le flux
    // d'energie / flux de volume Celles de quantité de mouvement avec une
    // valeur moyenne de u calculée par le flux de vitesse / flux de volume
    Real somme_flux_masse = 0.;
    Real somme_flux_volume = 0.;
    for (size_t imat = 0; imat < nbmat; imat++) somme_flux_volume += Flux[imat];

    if (std::fabs(somme_flux_volume) > flux_threshold ) {
      for (size_t imat = 0; imat < nbmat; imat++) {
        Flux[2. * nbmat + imat] =
            (Flux[2. * nbmat + imat] / somme_flux_volume) * Flux[nbmat + imat];
        somme_flux_masse += Flux[nbmat + imat];
      }

      Flux[3 * nbmat] = (Flux[3 * nbmat] / somme_flux_volume) *
                        somme_flux_masse;  // flux de quantité de mouvement x
      Flux[3 * nbmat + 1] =
          (Flux[3 * nbmat + 1] / somme_flux_volume) *
          somme_flux_masse;  // flux de quantité de mouvement y
      Flux[3 * nbmat + 2] =
          (Flux[3 * nbmat + 2] / somme_flux_volume) * somme_flux_masse;
      Flux[3 * nbmat + 3] =
          m_phi_lagrange[cell][3 * nbmat + 3] * somme_flux_volume;  // flux pour la pseudo VNR
    } else {
      Flux.fill(0.);
    }
  } else if (projectionPenteBorneDebarFix == 2) {
    // les flux d'energie massique se
    // deduisent des flux de volumes * densite (dans phi)
    Real somme_flux_masse = 0.;
    Real somme_flux_volume = 0.;
    for (size_t imat = 0; imat < nbmat; imat++) {
      Flux[2 * nbmat + imat] =
          m_phi_lagrange[cell][2 * nbmat + imat] *
          Flux[nbmat + imat];  // flux de masse energy de imat
      somme_flux_masse += Flux[nbmat + imat];
      somme_flux_volume += Flux[imat];
    }
    Flux[3 * nbmat] =
        m_phi_lagrange[cell][3 * nbmat] * somme_flux_masse;  // flux de quantité de mouvement x
    Flux[3 * nbmat + 1] = m_phi_lagrange[cell][3 * nbmat + 1] *
                          somme_flux_masse;  // flux de quantité de mouvement y
    Flux[3 * nbmat + 2] =
        m_phi_lagrange[cell][3 * nbmat + 2] * somme_flux_masse;  // flux d'energie cinetique
    Flux[3 * nbmat + 3] =
        m_phi_lagrange[cell][3 * nbmat + 3] * somme_flux_volume;  // flux pour la pseudo VNR
  }

  if (partie_positive_v == 0.) Flux.fill(0.);
  if (partie_positive_dual_v == 0.) Flux_dual.fill(0.);
  return;
}
//  *******************************************************************************
//  * \file computeRemapFlux
//  * \brief Calcul final des flux
//  *  cas classique : phi_face * v * dt * l * (n.ex ou n.ey)
//  *  phi_face contient la valeur de phi reconstruite à l'ordre voulu
//  *  cas pente-borne : phi_face * l * (n.ex ou n.ey)
//  *  phi_face contient deja le flux integre
//  * \param  
//  * \return valeur de phi 
//  *******************************************************************************
//  */
Real RemapADIService::computeRemapFlux(
        Integer projectionOrder, Integer projectionAvecPlateauPente,
        Real face_normal_velocity, Real3 face_normal,
        Real face_length, Real phi_face,
        Real3 outer_face_normal, Real3 exy, Real deltat_n) {
    
  if (projectionAvecPlateauPente == 0) {
// cas projection ordre 3 ou 1 ou 2 sans plateau pente (flux calculé ici
// avec phi_face)
    return (math::dot(outer_face_normal, exy) * face_normal_velocity * face_length *
            deltat_n * phi_face);
  } else {
    // cas projection ordre 2 avec plateau pente (flux dans la variable
    // phi_face)
    return (math::dot(outer_face_normal, exy) * face_length * phi_face);
  }
}
