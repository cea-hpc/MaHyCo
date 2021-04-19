// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "../MahycoModule.h"

/**
 *******************************************************************************
 * \file computeDualUremap1()
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
 *******************************************************************************
 */
void MahycoModule::computeDualUremap(Integer idir, String name)  {
  info() << " Entree dans computeDualUremap() pour la direction " << idir;
  Real deltat = m_global_deltat();
  NodeDirectionMng ndm(m_cartesian_mesh->nodeDirection(idir));
  m_dual_grad_phi.fill(Real3::null());
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
//       if (
//         // (m_node_coord[inode].x > 0.5 ) && (m_node_coord[inode].x < 0.51 ) )
//         (inode.localId() == 1809 )
//          || inode.localId() == 1810) 
//     {
//         info() << inode.localId() << " coord " << m_node_coord[inode] << " u_duale debut phase  " << idir << "   " << m_u_dual_lagrange[inode];
//     }
//       if (idir ==1) {
//       info() << backbacknode.localId() << " coord " << m_node_coord[backbacknode];
//       info() << backnode.localId() << " coord " << m_node_coord[backnode];
//       info() << inode.localId() << " coord " << m_node_coord[inode];
//       info() << frontnode.localId() << " coord " << m_node_coord[frontnode];
//       info() << frontfrontnode.localId() << " coord " << m_node_coord[frontfrontnode];
//       }
    }
    m_dual_grad_phi.synchronize();
  }
  String ajout_interne = "_INTERNE";
  name = name + ajout_interne;
  FaceGroup inner_dir_faces = mesh()->faceFamily()->findGroup(name);
  CellDirectionMng cdm(m_cartesian_mesh->cellDirection(idir));
  Integer nb_total_env = mm->environments().size();
  //m_back_flux_mass_env.fill(0.);
  ENUMERATE_NODE(inode, allNodes()){
    for (Integer index_env=0; index_env < nb_total_env; index_env++) {
      m_back_flux_mass_env[inode][index_env] =0.;
      m_front_flux_mass_env[inode][index_env] = 0.;
    }
    m_back_flux_mass[inode] = 0.;
    m_front_flux_mass[inode] = 0.;
  }
  // m_back_flux_mass.fill(0.);
  ENUMERATE_FACE(iface, inner_dir_faces) {
    Face face = *iface;       
    Cell cellb = face.backCell();
    Cell cellf = face.frontCell();
    Integer indexfacecellb(-1), indexfacecellf(-1);
    ENUMERATE_FACE(jface, cellb.faces()) { 
      if (jface.localId() == iface.localId()) indexfacecellb = jface.index(); 
    }
    ENUMERATE_FACE(jface, cellf.faces()) { 
      if (jface.localId() == iface.localId()) indexfacecellf = jface.index(); 
    }
    if ((indexfacecellb>6) || (indexfacecellf>6)) exit(1);
    if (! options()->projectionPenteBorne) {
      ENUMERATE_NODE(inode, face.nodes()) {
        // 4 faces dans une direction pour les noeuds --> 0.25  
        for (Integer index_env=0; index_env < nb_total_env; index_env++) { 
          m_back_flux_mass_env[inode][index_env] += 0.25 * m_flux_masse_face[cellb][indexfacecellb][index_env] * m_outer_face_normal[cellb][indexfacecellb][idir];
          m_front_flux_mass_env[inode][index_env] += 0.25 * m_flux_masse_face[cellf][indexfacecellf][index_env] * m_outer_face_normal[cellf][indexfacecellf][idir];
          // pour le flux total
          m_back_flux_mass[inode]  +=  m_back_flux_mass_env[inode][index_env];
          m_front_flux_mass[inode] +=  m_front_flux_mass_env[inode][index_env]; 
        }
//         if ( inode.localId() == 1809 || inode.localId() == 1810 ) {
//           info() <<  " bilan NODE " << inode.localId();
//           info() <<  " bilan " << m_back_flux_mass[inode] << " et " << m_front_flux_mass[inode];
//           info() <<  " 0env "  << m_back_flux_mass_env[inode][0] << " et no PB" << m_front_flux_mass_env[inode][0];
//           info() <<  " 1env "  << m_back_flux_mass_env[inode][1] << " et no PB" << m_front_flux_mass_env[inode][1];
//           info() << " flux front " << m_flux_masse_face[cellf][iface.localId()][0] << " face " << face.localId() << " " << m_flux_masse_face[cellf][iface.localId()][1]; 
//           info() << " flux back " << m_flux_masse_face[cellb][iface.localId()][0] << " face " << face.localId() << " " << m_flux_masse_face[cellb][iface.localId()][1]; 
//           info() << " normal back " <<  m_outer_face_normal[cellb][iface.localId()][idir] << " front " << m_outer_face_normal[cellf][iface.localId()][idir];
//           info() << " ----------------------------------------------------------------------------------";
//         }
      }
    } else {
      ENUMERATE_NODE(inode, face.nodes()) {
        // 4 faces dans une direction pour les noeuds --> 0.25  
        // recuperation du flux dual de masse calcule par le pente borne
        for (Integer index_env=0; index_env < nb_total_env; index_env++) { 
          m_back_flux_mass_env[inode][index_env] += 0.25 * m_dual_phi_flux[cellb][nb_total_env+index_env];
          m_front_flux_mass_env[inode][index_env] += 0.25 * m_dual_phi_flux[cellf][nb_total_env+index_env];
          // pour le flux total
          m_back_flux_mass[inode]  +=  m_back_flux_mass_env[inode][index_env];
          m_front_flux_mass[inode] +=  m_front_flux_mass_env[inode][index_env];
//           if ((m_node_coord[inode].x > 0.5 ) && (m_node_coord[inode].x < 0.52 ) ) {
//             info() <<  m_back_flux_mass[inode] << " et PB " << m_front_flux_mass[inode];
//             info() <<  m_back_flux_mass_env[inode][0] << " et PB" << m_front_flux_mass_env[inode][0];
//             info() <<  m_back_flux_mass_env[inode][1] << " et PB" << m_front_flux_mass_env[inode][1];
//           }
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
    // info() << " Passage avec le noeud " << node.localId() << " backnode " << backnode.localId() << " frontnode " << frontnode.localId();
    // info() << " u_duale av " << m_u_dual_lagrange[inode];
//     if ( inode.localId() == 1809 || inode.localId() == 1810 ) {
//           info() <<  " bilan NODE " << inode.localId();
//       info() << inode.localId() << " u_duale vitx AV " << m_u_dual_lagrange[inode][0] << " et m*vit " << m_node_mass[inode] * m_velocity[inode][0] ;
//       info() << inode.localId() << " u_duale vity AV " << m_u_dual_lagrange[inode][1] << " et m*vit " << m_node_mass[inode] * m_velocity[inode][1] ;
//       info() << inode.localId() << " u_duale vitz AV " << m_u_dual_lagrange[inode][2] << " et m*vit " << m_node_mass[inode] * m_velocity[inode][2] ;
//       info() << inode.localId() << " u_duale masse AV " << m_u_dual_lagrange[inode][3] << " et mass " << m_node_mass[inode];
//       info() <<  " bilan " << m_back_flux_mass[inode] << " et " << m_front_flux_mass[inode];
//           info() << " ----------------------------------------------------------------------------------";
//     }
    // flux de masse
    m_u_dual_lagrange[inode][3] += m_back_flux_mass[inode] - m_front_flux_mass[inode];
    
//     if ( inode.localId() == 1809 || inode.localId() == 1810 ) {
//       info() <<  " bilan NODE " << inode.localId();
//       info() << " m_front_flux_mass[inode] " << m_front_flux_mass[inode] << " et back " << m_back_flux_mass[inode] ;
//     }
    
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
    // info() << " puis frontupwind " <<  FrontupwindVelocity << " et back " << BackupwindVelocity;
//     if ((m_node_coord[inode].x > 0.5 ) && (m_node_coord[inode].x < 0.52 ) ) {
//         info() << BackupwindVelocity << " et " << FrontupwindVelocity;
//         info() <<  m_back_flux_mass[inode] << " et " << m_front_flux_mass[inode];
//         info() << " phiB " << m_phi_dual_lagrange[backnode];
//         info() << " phi " << m_phi_dual_lagrange[node];
//         info() << " phiF " << m_phi_dual_lagrange[frontnode];
//         info() << " grad phi " << m_dual_grad_phi[node];
//         info() << " uback " << uback;
//         info() << " ufront " << ufront;
//         info() << " Bnode " << backnode.localId() << " = " << m_node_coord[backnode][idir];
//         info() << " node " << node.localId() << " = " << m_node_coord[node][idir];
//         info() << " Fnode " << frontnode.localId() << " = " << m_node_coord[frontnode][idir];
//     }
    m_u_dual_lagrange[inode][0] += m_back_flux_mass[inode] * BackupwindVelocity[0] - 
        m_front_flux_mass[inode] * FrontupwindVelocity[0];
    m_u_dual_lagrange[inode][1] += m_back_flux_mass[inode] * BackupwindVelocity[1] - 
        m_front_flux_mass[inode] * FrontupwindVelocity[1];
    m_u_dual_lagrange[inode][2] += m_back_flux_mass[inode] * BackupwindVelocity[2] - 
        m_front_flux_mass[inode] * FrontupwindVelocity[2];
    // energie cinetique
    m_u_dual_lagrange[inode][4] += m_back_flux_mass[inode] * BackupwindEcin - 
        m_front_flux_mass[inode] * FrontupwindEcin;
        
//     if ( inode.localId() == 1809 || inode.localId() == 1810 ) {
//       info() <<  " bilan NODE " << inode.localId();
//       info() << inode.localId() << " u_duale vitx AP " << m_u_dual_lagrange[inode][0] ;
//       info() << inode.localId() << " u_duale vity AP " << m_u_dual_lagrange[inode][1] ;
//       info() << inode.localId() << " u_duale vitz AP " << m_u_dual_lagrange[inode][2] ;
//       info() << inode.localId() << " u_duale masse AP " << m_u_dual_lagrange[inode][3];
//       info() << inode.localId() << "BackupwindVelocity[1] " << BackupwindVelocity[1] << " " << m_back_flux_mass[inode];
//       info() << inode.localId() << " FrontupwindVelocity[1] " << FrontupwindVelocity[1] << " " << m_front_flux_mass[inode];
//       info() << " ----------------------------------------------------------------------------------";
//     }
        
    // recalcul des phi apres cette projection
    // phi vitesse
    m_phi_dual_lagrange[node][0] = m_u_dual_lagrange[inode][0] /  m_u_dual_lagrange[inode][3];
    m_phi_dual_lagrange[node][1] = m_u_dual_lagrange[inode][1] /  m_u_dual_lagrange[inode][3];
    m_phi_dual_lagrange[node][2] = m_u_dual_lagrange[inode][2] /  m_u_dual_lagrange[inode][3];
    // Phi masse nodale
    m_phi_dual_lagrange[node][3] = m_u_dual_lagrange[inode][3] ;
    // Phi energie
    m_phi_dual_lagrange[node][4] = m_u_dual_lagrange[inode][4] /  m_u_dual_lagrange[inode][3];
//     if (
//         // (m_node_coord[inode].x > 0.5 ) && (m_node_coord[inode].x < 0.51 ) )
//         (inode.localId() == 1809 )
//          || inode.localId() == 1810) 
//     {
//         info() << inode.localId() << " coord " << m_node_coord[inode] << " u_duale fin phase  " << idir << "   " << m_u_dual_lagrange[inode];
//     }
  }
  info() << " fin de la projection duale pour la direction " << idir;
}
/**
 *******************************************************************************
 * \file synchronizeUremap()
 * \brief phase de synchronisation des variables de projection dual apres projection
 * \return m_phi_dual_lagrange, m_u_dual_lagrange synchonise sur les mailles fantomes
 *******************************************************************************
 */
void MahycoModule::synchronizeDualUremap()  {
    debug() << " Entree dans synchronizeUremap()";
    m_phi_dual_lagrange.synchronize();
    m_u_dual_lagrange.synchronize();
}
