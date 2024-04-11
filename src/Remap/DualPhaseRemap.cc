// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapADIService.h"

/**
 *******************************************************************************
 * \file computeDualUremap()
 * \brief phase de projection duale
 *   Calcul du gradient aux noeuds 
 *   Calcul des flux de masse par environnement 
 *   m_back_flux_mass_env,  m_front_flux_mass_env
 *   et en moyenne
 *   m_back_flux_mass,  m_front_flux_mass
 *
 * \param idir : direction de la projection
 * \param name : nom du groupe de face suivant la direction demandée
 * \return m_dual_grad_phi
 * 
 * Résumé obtenu par ChatGPT
 * --------------------------
 * Cette fonction, nommée computeDualUremap(), est une méthode de la classe RemapADIService. 
 * Elle calcule les gradients aux nœuds ainsi que les flux de masse pour une projection duale. 
 * Plus précisément, elle calcule les flux de masse par environnement, m_back_flux_mass_env et m_front_flux_mass_env, 
 * ainsi que les flux de masse moyens m_back_flux_mass et m_front_flux_mass.
 * Les paramètres de cette fonction sont :
 * 
 * idir : la direction de la projection ;
 * nb_env : le nombre d'environnements considérés dans le calcul.
 * Le corps de la fonction commence par initialiser un certain nombre de variables et de tableaux, 
 * puis par calculer le gradient aux nœuds avec la méthode computeDualGradPhi(), si l'ordre de projection est supérieur à 1. 
 * Ensuite, la fonction calcule les flux de masse en parcourant toutes les faces de la direction considérée.
 * Les flux de masse sont calculés pour chaque nœud en sommant les flux associés aux cellules voisines. 
 * Finalement, les résultats sont synchronisés entre les processus.
 *******************************************************************************
 */
void RemapADIService::computeDualUremap(Integer idir, Integer nb_env)  {
    
  debug() << " Entree dans computeDualUremap() pour la direction " << idir;
  Real deltat = m_global_deltat();
  NodeDirectionMng ndm(m_cartesian_mesh->nodeDirection(idir));
  Real3 dirproj = {0.5 * (1-idir) * (2-idir), 
                   1.0 * idir * (2 -idir), 
                   -0.5 * idir * (1 - idir)};  
  m_dual_grad_phi.fill(0.0);
  if (options()->ordreProjection > 1) {
     ENUMERATE_NODE(inode, ndm.innerNodes()) {
      Node node = *inode;
      DirNode dir_node(ndm[inode]);
      Node backnode = dir_node.previous();
      DirNode dir_backnode(ndm[backnode]);
      Node backbacknode = dir_backnode.previous();
      Node frontnode = dir_node.next();
      DirNode dir_frontnode(ndm[frontnode]);
      Node frontfrontnode = dir_frontnode.next();

      computeDualGradPhi(node, frontfrontnode, frontnode, backnode, backbacknode, idir);   
    }
    m_dual_grad_phi.synchronize();
  }
  CellDirectionMng cdm(m_cartesian_mesh->cellDirection(idir));
  FaceDirectionMng fdm(m_cartesian_mesh->faceDirection(idir));
  //m_back_flux_mass_env.fill(0.);
  ENUMERATE_NODE(inode, allNodes()){
    for (Integer index_env=0; index_env < nb_env; index_env++) {
      m_back_flux_mass_env[inode][index_env] =0.;
      m_front_flux_mass_env[inode][index_env] = 0.;
    }
    m_back_flux_mass[inode] = 0.;
    m_front_flux_mass[inode] = 0.;
  }
  // 2 cellules dans une direction pour les noeuds --> 0.5  en 2D
  // 4 cellules dans une direction pour les noeuds --> 0.25  en 3D
  Real OneOverNbcell = ( mesh()->dimension() == 2 ? .5  : .25) ;
  ENUMERATE_FACE(iface, fdm.allFaces()) {
    Face face = *iface;       
    DirFace dir_face = fdm[face];
    Integer indexfacecellb(-1), indexfacecellf(-1);
    Cell cellb = dir_face.previousCell();
    if (cellb.localId() != -1) {
      ENUMERATE_NODE(inode, face.nodes()) {
        Node node = *inode;
        for (Integer index_env=0; index_env < nb_env; index_env++) { 
          // recuperation du flux dual de masse calcule par le pente borne ou 
          // dans le cas classique comme la somme des  flux à la cellule dans la direction donnée
          m_back_flux_mass_env[inode][index_env] +=   
          OneOverNbcell * m_dual_phi_flux[cellb][nb_env+index_env];
          // pour le flux total
          m_back_flux_mass[inode]  +=  OneOverNbcell * m_dual_phi_flux[cellb][nb_env+index_env];
         }
      }
    } 
    Cell cellf = dir_face.nextCell();    
    if (cellf.localId() != -1) {
      ENUMERATE_NODE(inode, face.nodes()) {
        Node node = *inode;
        for (Integer index_env=0; index_env < nb_env; index_env++) { 
          // recuperation du flux dual de masse calcule par le pente borne ou 
          // dans le cas classique comme la somme des  flux à la cellule dans la direction donnée
          m_front_flux_mass_env[inode][index_env] += 
          OneOverNbcell * m_dual_phi_flux[cellf][nb_env+index_env]; 
          // pour le flux total
          m_front_flux_mass[inode] += OneOverNbcell * m_dual_phi_flux[cellf][nb_env+index_env]; 
        }
      }
    }
  }
  
  // doit etre inutile si on augmente le nombre de mailles fantomes
  m_back_flux_mass_env.synchronize(); 
  m_front_flux_mass_env.synchronize();
  m_back_flux_mass.synchronize(); 
  m_front_flux_mass.synchronize();
    
  Real3 FrontupwindVelocity;
  Real3 BackupwindVelocity;
  Real FrontupwindEcin;
  Real BackupwindEcin;
  Integer order2 = options()->ordreProjection - 1;
  ENUMERATE_NODE(inode, ndm.innerNodes()) {
    Node node = *inode;
    DirNode dir_node(ndm[inode]);
    // info() << " Passage avec le noeud " << node.localId();
    Node backnode = dir_node.previous();
    Node frontnode = dir_node.next();
 
    // flux de masse
    m_u_dual_lagrange[inode][3] += m_back_flux_mass[inode] - m_front_flux_mass[inode];

    
    FrontupwindVelocity = 0.;
    BackupwindVelocity = 0.;
    // vitesse = vitesse(pNode) si FrontFluxMasse(pNode) > 0 et vitesse(voisin devant) sinon 
    Integer signfront;
    Real ufront = 0.5 * (m_phi_dual_lagrange[node][idir] + m_phi_dual_lagrange[frontnode][idir]);

    if (ufront > 0)  signfront = 1;
    else signfront = -1;
    if (m_front_flux_mass[inode] < 0.) {
      // vittese front upwind
      FrontupwindVelocity[0] = m_phi_dual_lagrange[frontnode][0]
      + order2 * ( 0.5 * m_dual_grad_phi[frontnode][0]  * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
      FrontupwindVelocity[1] = m_phi_dual_lagrange[frontnode][1]
      + order2 * ( 0.5 * m_dual_grad_phi[frontnode][1]  * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
      FrontupwindVelocity[2] = m_phi_dual_lagrange[frontnode][2]
      + order2 * ( 0.5 * m_dual_grad_phi[frontnode][2]  * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
      // energie cinetique
      FrontupwindEcin = m_phi_dual_lagrange[frontnode][4]
      + order2 * ( 0.5 * m_dual_grad_phi[frontnode][4] * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
    
    } else {
      // vittese front upwind
      FrontupwindVelocity[0] = m_phi_dual_lagrange[node][0]
      + order2 * ( 0.5 * m_dual_grad_phi[node][0]  * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
      FrontupwindVelocity[1] = m_phi_dual_lagrange[node][1]
      + order2 * ( 0.5 * m_dual_grad_phi[node][1]  * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
      FrontupwindVelocity[2] = m_phi_dual_lagrange[node][2]
      + order2 * ( 0.5 * m_dual_grad_phi[node][2]  * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
      // energie cinetique
      FrontupwindEcin = m_phi_dual_lagrange[node][4]
      + order2 * ( 0.5 * m_dual_grad_phi[node][4] * 
      (signfront * (m_node_coord[frontnode][idir] - m_node_coord[node][idir]) - deltat * ufront));
      
    }
    BackupwindVelocity = 0.;
    BackupwindEcin = 0.;
    // Backvitesse = vitesse(voisin amont)) si backFluxMasse(pNode) > 0 et vitesse(pNode) sinon
    Integer signback;
    Real  uback= 0.5 * (m_phi_dual_lagrange[backnode][idir] + m_phi_dual_lagrange[node][idir]);

    if (uback > 0) signback = 1;
    else signback= -1;

    if (m_back_flux_mass[inode] > 0.) {
      // vittese back upwind
      BackupwindVelocity[0] = m_phi_dual_lagrange[backnode][0]
      + order2 * ( 0.5 * m_dual_grad_phi[backnode][0]  * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
      BackupwindVelocity[1] = m_phi_dual_lagrange[backnode][1]
      + order2 * ( 0.5 * m_dual_grad_phi[backnode][1]  * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
      BackupwindVelocity[2] = m_phi_dual_lagrange[backnode][2]
      + order2 * ( 0.5 * m_dual_grad_phi[backnode][2]  * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
      // energie cinetique
      BackupwindEcin = m_phi_dual_lagrange[backnode][4]
      + order2 * ( 0.5 * m_dual_grad_phi[backnode][4] * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
      
    } else {
      // vittese back upwind
      BackupwindVelocity[0] = m_phi_dual_lagrange[node][0]
      + order2 * ( 0.5 * m_dual_grad_phi[node][0]  * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
      BackupwindVelocity[1] = m_phi_dual_lagrange[node][1]
      + order2 * ( 0.5 * m_dual_grad_phi[node][1]  * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
      BackupwindVelocity[2] = m_phi_dual_lagrange[node][2]
      + order2 * ( 0.5 * m_dual_grad_phi[node][2]  * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
      // energie cinetique
      BackupwindEcin = m_phi_dual_lagrange[node][4]
      + order2 * ( 0.5 * m_dual_grad_phi[node][4] * 
      (signback * (m_node_coord[node][idir] - m_node_coord[backnode][idir]) - deltat * uback));
      
    }
     
    m_u_dual_lagrange[inode][0] += m_back_flux_mass[inode] * BackupwindVelocity[0] - 
        m_front_flux_mass[inode] * FrontupwindVelocity[0];
    m_u_dual_lagrange[inode][1] += m_back_flux_mass[inode] * BackupwindVelocity[1] - 
        m_front_flux_mass[inode] * FrontupwindVelocity[1];
    m_u_dual_lagrange[inode][2] += m_back_flux_mass[inode] * BackupwindVelocity[2] - 
        m_front_flux_mass[inode] * FrontupwindVelocity[2];
    // energie cinetique
    m_u_dual_lagrange[inode][4] += m_back_flux_mass[inode] * BackupwindEcin - 
        m_front_flux_mass[inode] * FrontupwindEcin;
  }
  
  ENUMERATE_NODE(inode, ndm.innerNodes()) {
    Node node = *inode;
        
    // filtre des valeurs abherentes
    if (abs(m_u_dual_lagrange[node][0]) < m_arithmetic_thresold) m_u_dual_lagrange[node][0]=0.;
    if (abs(m_u_dual_lagrange[node][1]) < m_arithmetic_thresold) m_u_dual_lagrange[node][1]=0.;
    if (abs(m_u_dual_lagrange[node][2]) < m_arithmetic_thresold) m_u_dual_lagrange[node][2]=0.;
    // recalcul des phi apres cette projection
    // phi vitesse
    m_phi_dual_lagrange[node][0] = m_u_dual_lagrange[inode][0] /  m_u_dual_lagrange[inode][3];
    m_phi_dual_lagrange[node][1] = m_u_dual_lagrange[inode][1] /  m_u_dual_lagrange[inode][3];
    m_phi_dual_lagrange[node][2] = m_u_dual_lagrange[inode][2] /  m_u_dual_lagrange[inode][3];
    // Phi masse nodale
    m_phi_dual_lagrange[node][3] = m_u_dual_lagrange[inode][3] ;
    // Phi energie
    m_phi_dual_lagrange[node][4] = m_u_dual_lagrange[inode][4] /  m_u_dual_lagrange[inode][3];
  }
}
/*
 *******************************************************************************
 * \file synchronizeUremap()
 * \brief phase de synchronisation des variables de projection dual apres projection
 * \return m_phi_dual_lagrange, m_u_dual_lagrange synchonise sur les mailles fantomes
 *******************************************************************************
 */
void RemapADIService::synchronizeDualUremap()  {
    debug() << " Entree dans synchronizeUremap()";
    m_phi_dual_lagrange.synchronize();
    m_u_dual_lagrange.synchronize();
}
