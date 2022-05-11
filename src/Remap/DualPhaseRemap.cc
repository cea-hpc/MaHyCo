// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapADIService.h"
#include "accenv/AcceleratorUtils.h"
#include "cartesian/FactCartDirectionMng.h"


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
 *******************************************************************************
 */
void RemapADIService::computeDualUremap(Integer idir, Integer nb_env)  {
    
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans computeDualUremap() pour la direction " << idir;
  Real deltat = m_global_deltat();
  
  Real3 dirproj = {0.5 * (1-idir) * (2-idir), 
                   1.0 * idir * (2 -idir), 
                   -0.5 * idir * (1 - idir)};  
  
  auto queue = m_acc_env->newQueue();
  Cartesian::FactCartDirectionMng fact_cart(mesh());
  
  if (options()->ordreProjection > 1) {
    
    if (options()->projectionLimiteurId < minmodG) {
      
      // Il faut spécialiser pour les limiteurs < minmodG
      // Spécialisation
      switch (options()->projectionLimiteurId) {
        case minmod  : computeDualGradPhi_LimC<MinMod>   (idir); break;
        case superBee: computeDualGradPhi_LimC<SuperBee> (idir); break;
        case vanLeer : computeDualGradPhi_LimC<VanLeer>  (idir); break;
        default      : computeDualGradPhi_LimC<DefaultO1>(idir);
      }
    }
    else { // (options()->projectionLimiteurId > minmodG) PAS SUR GPU
      m_dual_grad_phi.fill(0.0);
      NodeDirectionMng ndm(m_cartesian_mesh->nodeDirection(idir));
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
    }
    m_dual_grad_phi.synchronize();
  }
  
#if 0
  ENUMERATE_NODE(inode, allNodes()){
    for (Integer index_env=0; index_env < nb_env; index_env++) {
      m_back_flux_mass_env[inode][index_env] =0.;
      m_front_flux_mass_env[inode][index_env] = 0.;
    }
    m_back_flux_mass[inode] = 0.;
    m_front_flux_mass[inode] = 0.;
  }
  
  // 2 cellules dans une direction pour les noeuds ==> 0.5 en 2D
  // 4 cellules dans une direction pour les noeuds ==> 0.25 en 3D
  Real oneovernbcell = ( mesh()->dimension() == 2 ? 0.5 : 0.25 );
  
  FaceDirectionMng fdm(m_cartesian_mesh->faceDirection(idir));
  
  ENUMERATE_FACE(iface, fdm.allFaces()) {
    Face face = *iface;       
    DirFace dir_face = fdm[face];
    Integer indexfacecellb(-1), indexfacecellf(-1);
    Cell cellb = dir_face.previousCell();
    if (cellb.localId() != -1) {
     ENUMERATE_FACE(jface, cellb.faces()) { 
      if (jface.localId() == iface.localId()) indexfacecellb = jface.index(); 
     }
     Real3 outer_face_normalb(m_outer_face_normal[cellb][indexfacecellb]);
//      Real outer_face_normal_dirb = math::dot(outer_face_normalb, dirproj);
     Real outer_face_normal_dirb = outer_face_normalb[idir];
     ENUMERATE_NODE(inode, face.nodes()) {
        for (Integer index_env=0; index_env < nb_env; index_env++) { 
        // recuperation du flux dual de masse calcule par le pente borne ou 
        // dans le cas classique comme 1 sur nb_cell de la somme des flux à la cellule dans la direction donnée,
        m_back_flux_mass_env[inode][index_env] += 
        oneovernbcell * m_dual_phi_flux[cellb][nb_env+index_env];
        // pour le flux total
        m_back_flux_mass[inode]  +=  oneovernbcell * m_dual_phi_flux[cellb][nb_env+index_env];
        }
     }
    }
    Cell cellf = dir_face.nextCell();
    if (cellf.localId() != -1) {
      ENUMERATE_FACE(jface, cellf.faces()) { 
       if (jface.localId() == iface.localId()) indexfacecellf = jface.index(); 
      }
      Real3 outer_face_normalf(m_outer_face_normal[cellf][indexfacecellf]);
//       Real outer_face_normal_dirf = math::dot(outer_face_normalf, dirproj);
      Real outer_face_normal_dirf = outer_face_normalf[idir];
      ENUMERATE_NODE(inode, face.nodes()) {
        for (Integer index_env=0; index_env < nb_env; index_env++) { 
        // recuperation du flux dual de masse calcule par le pente borne ou 
        // dans le cas classique comme comme 1 sur nb_cell de la somme des flux à la cellule dans la direction donnée,
        m_front_flux_mass_env[inode][index_env] += 
        oneovernbcell * m_dual_phi_flux[cellf][nb_env+index_env];
        // pour le flux total
        m_front_flux_mass[inode] += oneovernbcell * m_dual_phi_flux[cellf][nb_env+index_env];
        }
      }
    }
  }
  
#else 
   
  auto cfc = m_acc_env->connectivityView().cellFace();
  {
    auto command_f = makeCommand(queue);
    
    auto cart_fdm = fact_cart.faceDirection(idir);
    auto f2cid_stm = cart_fdm.face2CellIdStencil();
    auto face_group = cart_fdm.allFaces();
    
    auto in_dual_phi_flux      = ax::viewIn(command_f, m_dual_phi_flux    );
    auto in_outer_face_normal  = ax::viewIn(command_f, m_outer_face_normal);
    
    auto out_back_flux_contrib_env   = ax::viewOut(command_f, m_back_flux_contrib_env );
    auto out_front_flux_contrib_env  = ax::viewOut(command_f, m_front_flux_contrib_env );
    
    // 2 cellules dans une direction pour les noeuds ==> 0.5 en 2D
    // 4 cellules dans une direction pour les noeuds ==> 0.25 en 3D
    Real oneovernbcell = ( mesh()->dimension() == 2 ? 0.5 : 0.25 );
    
    command_f.addKernelName("fcontrib") << RUNCOMMAND_LOOP(iter, face_group.loopRanges()) {
      auto [fid, idx] = f2cid_stm.idIdx(iter); // id face + (i,j,k) face
      
      // Acces mailles gauche/droite 
      auto f2cid = f2cid_stm.face(fid, idx);
      CellLocalId backCid(f2cid.previousCell());
      CellLocalId frontCid(f2cid.nextCell());
      
      // Si face au bord gauche, on ne prend pas en compte la backCell
      if (!ItemId::null(backCid)) {
        Int16 index_face_backCid = 0;
        for( FaceLocalId backCid_fid : cfc.faces(backCid) ){
          if (backCid_fid == fid) {
            break;
          }
          index_face_backCid++;
        }
        Real3 outer_face_normalb(in_outer_face_normal[backCid][index_face_backCid]);
        Real outer_face_normal_dirb = outer_face_normalb[idir];

        for (Integer index_env=0; index_env < nb_env; index_env++) {
          out_back_flux_contrib_env[backCid][index_env] = oneovernbcell * in_dual_phi_flux[backCid][nb_env+index_env];
        }
      }
      
      // Si face au bord droit, on ne prend pas en compte la frontCell
      if (!ItemId::null(frontCid)) {
        Int16 index_face_frontCid = 0;
        for( FaceLocalId frontCid_fid : cfc.faces(frontCid) ){
          if (frontCid_fid == fid) {
            break;
          }
          index_face_frontCid++;
        }
        Real3 outer_face_normalf(in_outer_face_normal[frontCid][index_face_frontCid]);
        Real outer_face_normal_dirf = outer_face_normalf[idir];
        
        for (Integer index_env=0; index_env < nb_env; index_env++) {
          out_front_flux_contrib_env[frontCid][index_env] = oneovernbcell * in_dual_phi_flux[frontCid][nb_env+index_env];
        }
      }
    };
  }
  
  if (mesh()->dimension() == 2)
  {
    auto command = makeCommand(queue);
    
    auto in_back_flux_contrib_env   = ax::viewIn(command, m_back_flux_contrib_env );
    auto in_front_flux_contrib_env  = ax::viewIn(command, m_front_flux_contrib_env );
 
    auto inout_back_flux_mass      = ax::viewInOut(command,m_back_flux_mass     );
    auto inout_front_flux_mass     = ax::viewInOut(command,m_front_flux_mass    );
    auto inout_back_flux_mass_env  = ax::viewInOut(command,m_back_flux_mass_env );
    auto inout_front_flux_mass_env = ax::viewInOut(command,m_front_flux_mass_env);
    
    auto cart_cell_dm = fact_cart.cellDirection(idir);
    auto c2cid_stm = cart_cell_dm.cell2CellIdStencil();

    auto node_dm = fact_cart.nodeDirection(idir);
    auto n2nid_stm = node_dm.node2NodeIdStencil();
    auto node_group = node_dm.allNodes();

    // Direction transverse en 2D
    Integer dir1 = (idir+1)%2;
    Integer ncells0 = fact_cart.cartesianGrid()->cartNumCell().nbItemDir(0);
    Integer ncells1 = fact_cart.cartesianGrid()->cartNumCell().nbItemDir(1);

    command.addKernelName("ncontrib") << RUNCOMMAND_LOOP(iter, node_group.loopRanges()) {
      auto [nid, idx] = n2nid_stm.idIdx(iter); // node id + (i,j,k) du noeud

      for (Integer index_env=0; index_env < nb_env; index_env++) {
        inout_back_flux_mass_env[nid][index_env] =0.;
        inout_front_flux_mass_env[nid][index_env] = 0.;
      }
      inout_back_flux_mass[nid] = 0.;
      inout_front_flux_mass[nid] = 0.;

      // 2 Mailles "back" dans la direction en 2D
      for(Integer sht1=-1 ; sht1<=0 ; ++sht1) {
        IdxType cellb_idx=idx;
        cellb_idx[idir]-=1; // back
        cellb_idx[dir1]+=sht1;

        // La maille existe-t-elle ?
        if ((cellb_idx[0]>=0 && cellb_idx[0]<ncells0) &&
            (cellb_idx[1]>=0 && cellb_idx[1]<ncells1)) {
          CellLocalId backCid{c2cid_stm.id(cellb_idx[0],cellb_idx[1],cellb_idx[2])};

          for (Integer index_env=0; index_env < nb_env; index_env++) { 
            inout_back_flux_mass_env[nid][index_env] += in_back_flux_contrib_env[backCid][index_env];
            inout_back_flux_mass    [nid]            += in_back_flux_contrib_env[backCid][index_env];
          }
        }
      }

      // 2 Mailles "front" dans la direction en 2D
      for(Integer sht1=-1 ; sht1<=0 ; ++sht1) {
        IdxType cellf_idx=idx;
        // cellf_idx[idir]+=0; // front
        cellf_idx[dir1]+=sht1;

        // La maille existe-t-elle ?
        if ((cellf_idx[0]>=0 && cellf_idx[0]<ncells0) &&
            (cellf_idx[1]>=0 && cellf_idx[1]<ncells1)) {
          CellLocalId frontCid{c2cid_stm.id(cellf_idx[0],cellf_idx[1],cellf_idx[2])};

          for (Integer index_env=0; index_env < nb_env; index_env++) { 
            inout_front_flux_mass_env[nid][index_env] += in_front_flux_contrib_env[frontCid][index_env];
            inout_front_flux_mass    [nid]            += in_front_flux_contrib_env[frontCid][index_env];
          }
        }
      }
    };
  }
  else if (mesh()->dimension() == 3)
  {
    auto command = makeCommand(queue);
    
    auto in_back_flux_contrib_env   = ax::viewIn(command, m_back_flux_contrib_env );
    auto in_front_flux_contrib_env  = ax::viewIn(command, m_front_flux_contrib_env );
 
    auto inout_back_flux_mass      = ax::viewInOut(command,m_back_flux_mass     );
    auto inout_front_flux_mass     = ax::viewInOut(command,m_front_flux_mass    );
    auto inout_back_flux_mass_env  = ax::viewInOut(command,m_back_flux_mass_env );
    auto inout_front_flux_mass_env = ax::viewInOut(command,m_front_flux_mass_env);
    
    auto cart_cell_dm = fact_cart.cellDirection(idir);
    auto c2cid_stm = cart_cell_dm.cell2CellIdStencil();

    auto node_dm = fact_cart.nodeDirection(idir);
    auto n2nid_stm = node_dm.node2NodeIdStencil();
    auto node_group = node_dm.allNodes();

    // Directions transverses en 3D
    Integer dir1 = (idir+1)%3;
    Integer dir2 = (idir+2)%3;
    Integer ncells0 = fact_cart.cartesianGrid()->cartNumCell().nbItemDir(0);
    Integer ncells1 = fact_cart.cartesianGrid()->cartNumCell().nbItemDir(1);
    Integer ncells2 = fact_cart.cartesianGrid()->cartNumCell().nbItemDir(2);

    command.addKernelName("ncontrib") << RUNCOMMAND_LOOP(iter, node_group.loopRanges()) {
      auto [nid, idx] = n2nid_stm.idIdx(iter); // node id + (i,j,k) du noeud

      for (Integer index_env=0; index_env < nb_env; index_env++) {
        inout_back_flux_mass_env[nid][index_env] =0.;
        inout_front_flux_mass_env[nid][index_env] = 0.;
      }
      inout_back_flux_mass[nid] = 0.;
      inout_front_flux_mass[nid] = 0.;

      // 4 Mailles "back" dans la direction en 3D
      for(Integer sht1=-1 ; sht1<=0 ; ++sht1) {
        for(Integer sht2=-1 ; sht2<=0 ; ++sht2) {
          IdxType cellb_idx=idx;
          cellb_idx[idir]-=1; // back
          cellb_idx[dir1]+=sht1;
          cellb_idx[dir2]+=sht2;

          // La maille existe-t-elle ?
          if ((cellb_idx[0]>=0 && cellb_idx[0]<ncells0) &&
              (cellb_idx[1]>=0 && cellb_idx[1]<ncells1) &&
              (cellb_idx[2]>=0 && cellb_idx[2]<ncells2)) {
            CellLocalId backCid{c2cid_stm.id(cellb_idx[0],cellb_idx[1],cellb_idx[2])};

            for (Integer index_env=0; index_env < nb_env; index_env++) { 
              inout_back_flux_mass_env[nid][index_env] += in_back_flux_contrib_env[backCid][index_env];
              inout_back_flux_mass    [nid]            += in_back_flux_contrib_env[backCid][index_env];
            }
          }
        }
      }

      // 4 Mailles "front" dans la direction en 3D
      for(Integer sht1=-1 ; sht1<=0 ; ++sht1) {
        for(Integer sht2=-1 ; sht2<=0 ; ++sht2) {
          IdxType cellf_idx=idx;
          // cellf_idx[idir]+=0; // front
          cellf_idx[dir1]+=sht1;
          cellf_idx[dir2]+=sht2;

          // La maille existe-t-elle ?
          if ((cellf_idx[0]>=0 && cellf_idx[0]<ncells0) &&
              (cellf_idx[1]>=0 && cellf_idx[1]<ncells1) &&
              (cellf_idx[2]>=0 && cellf_idx[2]<ncells2)) {
            CellLocalId frontCid{c2cid_stm.id(cellf_idx[0],cellf_idx[1],cellf_idx[2])};

            for (Integer index_env=0; index_env < nb_env; index_env++) { 
              inout_front_flux_mass_env[nid][index_env] += in_front_flux_contrib_env[frontCid][index_env];
              inout_front_flux_mass    [nid]            += in_front_flux_contrib_env[frontCid][index_env];
            }
          }
        }
      }
    };
  }
#endif

  // doit etre inutile si on augmente le nombre de mailles fantomes
  m_back_flux_mass_env.synchronize(); 
  m_front_flux_mass_env.synchronize();
  m_back_flux_mass.synchronize(); 
  m_front_flux_mass.synchronize();
    
#if 1
  {
    Integer order2 = options()->ordreProjection - 1;
    Real thresold = m_arithmetic_thresold;
    
    auto command = makeCommand(queue);
    
    auto cart_ndm = fact_cart.nodeDirection(idir);
    auto n2nid_stm = cart_ndm.node2NodeIdStencil();
    
    auto node_group = cart_ndm.innerNodes();
    
    auto in_back_flux_mass  = ax::viewIn(command, m_back_flux_mass );
    auto in_front_flux_mass = ax::viewIn(command, m_front_flux_mass);
    auto in_dual_grad_phi   = ax::viewIn(command, m_dual_grad_phi  );
    auto in_node_coord      = ax::viewIn(command, m_node_coord     );
    
    auto inout_u_dual_lagrange   = ax::viewInOut(command, m_u_dual_lagrange  );
    auto inout_phi_dual_lagrange = ax::viewInOut(command, m_phi_dual_lagrange);
    
    
    command << RUNCOMMAND_LOOP(iter, node_group.loopRanges()) {
      auto [nid, idx] = n2nid_stm.idIdx(iter); // id maille + (i,j,k) maille
      
      // Acces noeuds gauche/droite qui existent forcement
      auto n2nid = n2nid_stm.stencilNode<1>(nid, idx);
      
      NodeLocalId backNid(n2nid.previousId()); // back node
      NodeLocalId frontNid(n2nid.nextId()); // front node
      
      // flux de masse
      inout_u_dual_lagrange[nid][3] += in_back_flux_mass[nid] - in_front_flux_mass[nid];
      
      Real3 FrontupwindVelocity(0., 0., 0.);
      Real3 BackupwindVelocity (0., 0., 0.);
      Real  FrontupwindEcin = 0.;
      Real  BackupwindEcin  = 0.;
      
      // vitesse = vitesse(pNode) si FrontFluxMasse(pNode) > 0 et vitesse(voisin devant) sinon 
      Integer signfront;
      Real ufront = 0.5 * (inout_phi_dual_lagrange[nid][idir] + inout_phi_dual_lagrange[nid][idir]);

      if (ufront > 0)  signfront = 1;
      else signfront = -1;
      
      if (in_front_flux_mass[nid] < 0.) {
        // vittese front upwind
        FrontupwindVelocity[0] = inout_phi_dual_lagrange[frontNid][0]
        + order2 * ( 0.5 * in_dual_grad_phi[frontNid][0]  * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
        
        FrontupwindVelocity[1] = inout_phi_dual_lagrange[frontNid][1]
        + order2 * ( 0.5 * in_dual_grad_phi[frontNid][1]  * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
        
        FrontupwindVelocity[2] = inout_phi_dual_lagrange[frontNid][2]
        + order2 * ( 0.5 * in_dual_grad_phi[frontNid][2]  * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
        
        // energie cinetique
        FrontupwindEcin = inout_phi_dual_lagrange[frontNid][4]
        + order2 * ( 0.5 * in_dual_grad_phi[frontNid][4] * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
        
      } else {
        // vittese front upwind
        FrontupwindVelocity[0] = inout_phi_dual_lagrange[nid][0]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][0]  * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
        
        FrontupwindVelocity[1] = inout_phi_dual_lagrange[nid][1]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][1]  * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
        
        FrontupwindVelocity[2] = inout_phi_dual_lagrange[nid][2]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][2]  * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
        
        // energie cinetique
        FrontupwindEcin = inout_phi_dual_lagrange[nid][4]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][4] * 
        (signfront * (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]) - deltat * ufront));
      }
       
      // Backvitesse = vitesse(voisin amont)) si backFluxMasse(pNode) > 0 et vitesse(pNode) sinon
      Integer signback;
      Real  uback= 0.5 * (inout_phi_dual_lagrange[backNid][idir] + inout_phi_dual_lagrange[nid][idir]);
      
      if (uback > 0) signback = 1;
      else signback= -1;
      
      if (in_back_flux_mass[nid] > 0.) {
        // vittese back upwind
        BackupwindVelocity[0] = inout_phi_dual_lagrange[backNid][0]
        + order2 * ( 0.5 * in_dual_grad_phi[backNid][0]  * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
        
        BackupwindVelocity[1] = inout_phi_dual_lagrange[backNid][1]
        + order2 * ( 0.5 * in_dual_grad_phi[backNid][1]  * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
        
        BackupwindVelocity[2] = inout_phi_dual_lagrange[backNid][2]
        + order2 * ( 0.5 * in_dual_grad_phi[backNid][2]  * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
        
        // energie cinetique
        BackupwindEcin = inout_phi_dual_lagrange[backNid][4]
        + order2 * ( 0.5 * in_dual_grad_phi[backNid][4] * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
        
        
      } else {
        // vittese back upwind
        BackupwindVelocity[0] = inout_phi_dual_lagrange[nid][0]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][0]  * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
        
        BackupwindVelocity[1] = inout_phi_dual_lagrange[nid][1]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][1]  * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
        
        BackupwindVelocity[2] = inout_phi_dual_lagrange[nid][2]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][2]  * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
        
        // energie cinetique
        BackupwindEcin = inout_phi_dual_lagrange[nid][4]
        + order2 * ( 0.5 * in_dual_grad_phi[nid][4] * 
        (signback * (in_node_coord[nid][idir] - in_node_coord[backNid][idir]) - deltat * uback));
      }
      
      inout_u_dual_lagrange[nid][0] += in_back_flux_mass[nid] * BackupwindVelocity[0] - 
      in_front_flux_mass[nid] * FrontupwindVelocity[0];
      inout_u_dual_lagrange[nid][1] += in_back_flux_mass[nid] * BackupwindVelocity[1] - 
      in_front_flux_mass[nid] * FrontupwindVelocity[1];
      inout_u_dual_lagrange[nid][2] += in_back_flux_mass[nid] * BackupwindVelocity[2] - 
      in_front_flux_mass[nid] * FrontupwindVelocity[2];
      // energie cinetique
      inout_u_dual_lagrange[nid][4] += in_back_flux_mass[nid] * BackupwindEcin - 
      in_front_flux_mass[nid] * FrontupwindEcin;
    };
  }
  
  {
    auto command = makeCommand(queue);
    
    Real thresold = m_arithmetic_thresold;
    
    auto cart_ndm = fact_cart.nodeDirection(idir);
    auto n2nid_stm = cart_ndm.node2NodeIdStencil();
    
    auto node_group = cart_ndm.innerNodes();
    
    auto inout_u_dual_lagrange   = ax::viewInOut(command, m_u_dual_lagrange  );
    auto inout_phi_dual_lagrange = ax::viewInOut(command, m_phi_dual_lagrange);
    
    command << RUNCOMMAND_LOOP(iter, node_group.loopRanges()) {
      auto [nid, idx] = n2nid_stm.idIdx(iter); // id maille + (i,j,k) maille      
      // filtre des valeurs abherentes
      if (abs(inout_u_dual_lagrange[nid][0]) < thresold) inout_u_dual_lagrange[nid][0]=0.;
      if (abs(inout_u_dual_lagrange[nid][1]) < thresold) inout_u_dual_lagrange[nid][1]=0.;
      if (abs(inout_u_dual_lagrange[nid][2]) < thresold) inout_u_dual_lagrange[nid][2]=0.;
      // recalcul des phi apres cette projection
      // phi vitesse
      inout_phi_dual_lagrange[nid][0] = inout_u_dual_lagrange[nid][0] /  inout_u_dual_lagrange[nid][3];
      inout_phi_dual_lagrange[nid][1] = inout_u_dual_lagrange[nid][1] /  inout_u_dual_lagrange[nid][3];
      inout_phi_dual_lagrange[nid][2] = inout_u_dual_lagrange[nid][2] /  inout_u_dual_lagrange[nid][3];
      // Phi masse nodale
      inout_phi_dual_lagrange[nid][3] = inout_u_dual_lagrange[nid][3] ;
      // Phi energie
      inout_phi_dual_lagrange[nid][4] = inout_u_dual_lagrange[nid][4] /  inout_u_dual_lagrange[nid][3];
    };
  }
#else
  
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
  
  

#endif

  PROF_ACC_END;
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


/**
 * ******************************************************************************
 * \file computeDualGradPhi_LimC()
 * \brief Spécialisation de computeDualGradPhi
 *        et options()->projectionLimiteurId < minmodG (limiteur classique)
 * \param
 * \return 
 *******************************************************************************
 */
template<typename LimType>
void RemapADIService::
computeDualGradPhi_LimC(Integer idir) {
  PROF_ACC_BEGIN(__FUNCTION__);
  
  Cartesian::FactCartDirectionMng fact_cart(mesh());
  
  auto queue = m_acc_env->newQueue();
  {
    auto command = makeCommand(queue);
    
    auto cart_ndm = fact_cart.nodeDirection(idir);
    auto n2nid_stm = cart_ndm.node2NodeIdStencil();
    
    auto node_group = cart_ndm.innerNodes();
    
    auto in_phi_dual_lagrange = ax::viewIn(command, m_phi_dual_lagrange);
    auto in_node_coord        = ax::viewIn(command, m_node_coord);
    
    auto out_dual_grad_phi = ax::viewOut(command, m_dual_grad_phi);
    
    command << RUNCOMMAND_LOOP(iter, node_group.loopRanges()) {
      auto [nid, idx] = n2nid_stm.idIdx(iter); // id maille + (i,j,k) maille
      
      // Acces noeuds gauche/droite qui existent forcement
      auto n2nid = n2nid_stm.stencilNode<2>(nid, idx);
      
      NodeLocalId backNid(n2nid.previousId()); // back node
      NodeLocalId frontNid(n2nid.nextId()); // front node
      NodeLocalId backbackNid(n2nid.prev_previousId()); // back back node
      NodeLocalId frontfrontNid(n2nid.next_nextId()); // front front node
      
      Real3 grad_front(0. , 0. , 0.);
      Real3 grad_back (0. , 0. , 0.);
      
      // gradient vitesse X selon la direction idir
      grad_front[0] = (in_phi_dual_lagrange[frontNid][0] - in_phi_dual_lagrange[nid][0]) /
      (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]);
      
      grad_back[0] = (in_phi_dual_lagrange[nid][0] - in_phi_dual_lagrange[backNid][0]) /
      (in_node_coord[nid][idir] - in_node_coord[backNid][idir]);
      
      // gradient vitesse Y selon la direction idir
      grad_front[1] = (in_phi_dual_lagrange[frontNid][1] - in_phi_dual_lagrange[nid][1]) /
      (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]);
      
      grad_back[1] = (in_phi_dual_lagrange[nid][1] - in_phi_dual_lagrange[backNid][1]) /
      (in_node_coord[nid][idir] - in_node_coord[backNid][idir]);
      
      // gradient vitesse Z selon la direction idir
      grad_front[2] = (in_phi_dual_lagrange[frontNid][2] - in_phi_dual_lagrange[nid][2]) /
      (in_node_coord[frontNid][idir] - in_node_coord[nid][idir]);
      
      grad_back[2] = (in_phi_dual_lagrange[nid][2] - in_phi_dual_lagrange[backNid][2]) /
      (in_node_coord[nid][idir] - in_node_coord[backNid][idir]);
      
      // largeurs des mailles duales
      Real hmoins, h0, hplus;
      h0 = 0.5 * (in_node_coord[frontNid][idir]- in_node_coord[backNid][idir]);
      if (ItemId::null(backbackNid)) {
        hmoins = 0.;
        hplus = 0.5 * (in_node_coord[frontfrontNid][idir]- in_node_coord[nid][idir]);
      } else if (ItemId::null(frontfrontNid)) {
        hplus = 0.;
        hmoins = 0.5 * (in_node_coord[nid][idir] - in_node_coord[backbackNid][idir]);
      } else {
        hmoins = 0.5 * (in_node_coord[nid][idir] - in_node_coord[backbackNid][idir]);
        hplus = 0.5 * (in_node_coord[frontfrontNid][idir]- in_node_coord[nid][idir]);
      }
      
      // info() << " Passage gradient limite Classique ";
      for (Integer ivar = 0; ivar < 3; ivar++) {
        out_dual_grad_phi[nid][ivar] = 0.;
        if (grad_back[ivar] != 0. ) 
          out_dual_grad_phi[nid][ivar] += 0.5 * LimType::fluxLimiter( grad_front[ivar] / grad_back[ivar] ) * grad_back[ivar];
        if (grad_front[ivar] !=0.)
          out_dual_grad_phi[nid][ivar] += 0.5 * LimType::fluxLimiter( grad_back[ivar] / grad_front[ivar] ) * grad_front[ivar];
      }
    };
  }
  
  
  
//   NodeDirectionMng ndm(m_cartesian_mesh->nodeDirection(idir));
//   
//   ENUMERATE_NODE(inode, ndm.innerNodes()) {
//     
//     Node node = *inode;
//     DirNode dir_node(ndm[inode]);
//     Node backnode = dir_node.previous();
//     DirNode dir_backnode(ndm[backnode]);
//     Node backbacknode = dir_backnode.previous();
//     Node frontnode = dir_node.next();
//     DirNode dir_frontnode(ndm[frontnode]);
//     Node frontfrontnode = dir_frontnode.next();
//     
//     Real3 grad_front(0. , 0. , 0.);
//     Real3 grad_back(0. , 0. , 0.);
//     
//     
//     // gradient vitesse X selon la direction idir
//     grad_front[0] = (m_phi_dual_lagrange[frontnode][0] - m_phi_dual_lagrange[node][0]) /
//     (m_node_coord[frontnode][idir] - m_node_coord[node][idir]);
//     
//     grad_back[0] = (m_phi_dual_lagrange[node][0] - m_phi_dual_lagrange[backnode][0]) /
//     (m_node_coord[node][idir] - m_node_coord[backnode][idir]);
//     
//     // gradient vitesse Y selon la direction idir
//     grad_front[1] = (m_phi_dual_lagrange[frontnode][1] - m_phi_dual_lagrange[node][1]) /
//     (m_node_coord[frontnode][idir] - m_node_coord[node][idir]);
//     
//     grad_back[1] = (m_phi_dual_lagrange[node][1] - m_phi_dual_lagrange[backnode][1]) /
//     (m_node_coord[node][idir] - m_node_coord[backnode][idir]);
//     
//     // gradient vitesse Z selon la direction idir
//     grad_front[2] = (m_phi_dual_lagrange[frontnode][2] - m_phi_dual_lagrange[node][2]) /
//     (m_node_coord[frontnode][idir] - m_node_coord[node][idir]);
//     
//     grad_back[2] = (m_phi_dual_lagrange[node][2] - m_phi_dual_lagrange[backnode][2]) /
//     (m_node_coord[node][idir] - m_node_coord[backnode][idir]);
//     
//     // largeurs des mailles duales
//     Real hmoins, h0, hplus;
//     h0 = 0.5 * (m_node_coord[frontnode][idir]- m_node_coord[backnode][idir]);
//     if (backbacknode.localId() == -1) {
//       hmoins = 0.;
//       hplus = 0.5 * (m_node_coord[frontfrontnode][idir]- m_node_coord[node][idir]);
//     } else if (frontfrontnode.localId() == -1) {
//       hplus = 0.;
//       hmoins = 0.5 * (m_node_coord[node][idir] - m_node_coord[backbacknode][idir]);
//     } else {
//       hmoins = 0.5 * (m_node_coord[node][idir] - m_node_coord[backbacknode][idir]);
//       hplus = 0.5 * (m_node_coord[frontfrontnode][idir]- m_node_coord[node][idir]);
//     }
//     
//     // info() << " Passage gradient limite Classique ";
//     for (Integer ivar = 0; ivar < 3; ivar++) {
//                 m_dual_grad_phi[node][ivar] = 0.;
//                 if (grad_back[ivar] != 0. ) 
//                   m_dual_grad_phi[node][ivar] += 0.5 * LimType::fluxLimiter( grad_front[ivar] / grad_back[ivar] ) * grad_back[ivar];
//                 if (grad_front[ivar] !=0.)
//                   m_dual_grad_phi[node][ivar] += 0.5 * LimType::fluxLimiter( grad_back[ivar] / grad_front[ivar] ) * grad_front[ivar];
//     }
//   }
  PROF_ACC_END;
}
