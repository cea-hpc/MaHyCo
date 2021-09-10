// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MahycoModule::
_computeNodeIndexInCells() {
  debug() << "_computeNodeIndexInCells";
  // Un noeud est connecté au maximum à MAX_NODE_CELL mailles
  // Calcul pour chaque noeud son index dans chacune des
  // mailles à laquelle il est connecté.
  NodeGroup nodes = allNodes();
  Integer nb_node = nodes.size();
  m_node_index_in_cells.resize(MAX_NODE_CELL*nb_node);
  m_node_index_in_cells.fill(-1);
  auto node_cell_cty = m_connectivity_view.nodeCell();
  auto cell_node_cty = m_connectivity_view.cellNode();
  ENUMERATE_NODE(inode,nodes){
    NodeLocalId node = *inode;
    Int32 index = 0; 
    Int32 first_pos = node.localId() * MAX_NODE_CELL;
    for( CellLocalId cell : node_cell_cty.cells(node) ){
      Int16 node_index_in_cell = 0; 
      for( NodeLocalId cell_node : cell_node_cty.nodes(cell) ){
        if (cell_node==node)
          break;
        ++node_index_in_cell;
      }    
      m_node_index_in_cells[first_pos + index] = node_index_in_cell;
      ++index;
    }    
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
_computeNodeIndexInFaces()
{
  info() << "ComputeNodeIndexInFaces";
  // Un noeud est connecté au maximum à MAX_NODE_FACE facettes
  // Calcul pour chaque noeud son index dans chacune des
  // faces à laquelle il est connecté.
  NodeGroup nodes = allNodes();
  Integer nb_node = nodes.size();
  m_node_index_in_faces.resize(MAX_NODE_FACE*nb_node);
  m_node_index_in_faces.fill(-1);
  auto node_face_cty = m_connectivity_view.nodeFace();
  auto face_node_cty = m_connectivity_view.faceNode();
  ENUMERATE_NODE(inode,nodes){
    NodeLocalId node = *inode;
    Int32 index = 0;
    Int32 first_pos = node.localId() * MAX_NODE_FACE;
    for( FaceLocalId face : node_face_cty.faces(node) ){
      Int16 node_index_in_face = 0;
      for( NodeLocalId face_node : face_node_cty.nodes(face) ){
        if (face_node==node)
          break;
        ++node_index_in_face;
      }
      m_node_index_in_faces[first_pos + index] = node_index_in_face;
      ++index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/* A appeler par hydroStartInit et par hydroContinueInit pour préparer les   */
/* données pour les accélérateurs                                            */
/*---------------------------------------------------------------------------*/
void MahycoModule::
_initMeshForAcc() {
  debug() << "_initMeshForAcc";

  m_connectivity_view.setMesh(this->mesh());

  // Permet la lecture des cqs quand on boucle sur les noeuds
  _computeNodeIndexInCells();
  _computeNodeIndexInFaces();

#ifdef ARCANE_HAS_CUDA
  // "Conseils" utilisation de la mémoire unifiée
  int device = -1;
  cudaGetDevice(&device);

  mem_adv_set_read_mostly(m_node_index_in_cells.view(), device);
  mem_adv_set_read_mostly(m_node_index_in_faces.view(), device);
  
  // CellLocalId
  mem_adv_set_read_mostly(allCells().view().localIds(), device);
  mem_adv_set_read_mostly(ownCells().view().localIds(), device);

  // NodeLocalId
  mem_adv_set_read_mostly(allNodes().view().localIds(), device);
  mem_adv_set_read_mostly(ownNodes().view().localIds(), device);

  // FaceLocalId
  mem_adv_set_read_mostly(allFaces().view().localIds(), device);
  mem_adv_set_read_mostly(ownFaces().view().localIds(), device);

  for(Integer dir(0) ; dir < mesh()->dimension() ; ++dir) {
    auto cell_dm = m_cartesian_mesh->cellDirection(dir);

    // "Conseils" accés mémoire
    mem_adv_set_read_mostly(cell_dm.innerCells().view().localIds(), device);
    mem_adv_set_read_mostly(cell_dm.outerCells().view().localIds(), device);
    mem_adv_set_read_mostly(cell_dm.allCells().view().localIds(), device);
  }
#endif
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::hydroStartInitEnvAndMat()
{
  info() << " Preparation des env  ";
  IMeshMaterialMng* m_material_mng;
  mm = IMeshMaterialMng::getReference(defaultMesh());
  // lecture des environnements
  MeshBlockBuildInfo mbbi("BLOCK1",allCells());
  UniqueArray<IMeshEnvironment*> saved_envs;

  info() << "Lit les infos des matériaux du JDD ";
  // Lit les infos des matériaux du JDD et les enregistre dans le gestionnaire
  for( Integer i=0,n=options()->material().size(); i<n; ++i ){
    String mat_name = options()->material[i].name;
    info() << "Found material name=" << mat_name;
    mm->registerMaterialInfo(mat_name);
  }

  for( Integer i=0,n=options()->environment().size(); i<n; ++i ){
    String env_name = options()->environment[i].name;
    info() << "Found environment name=" << env_name;
    Materials::MeshEnvironmentBuildInfo env_build(env_name);
    info() << " nbmat " << options()->environment[i].material.size();
    for( Integer k=0,kn=options()->environment[i].material.size(); k<kn; ++k ){
      String mat_name = options()->environment[i].material[k];
      info() << " k = " << k << " Add material " << mat_name << " for environment " << env_name;
      env_build.addMaterial(mat_name);
    }
    info() << "Materiau cree";
    IMeshEnvironment* env = mm->createEnvironment(env_build);
    info() << "Environment cree";
    saved_envs.add(env);
    // Le bloc ne contient que 2 milieux
    if (i<2){
      info() << "Add environment " << env_name << " to block1";
      mbbi.addEnvironment(env);
    }
  }
  
  
  info() << " Rangement des mailles  ";
  
  IMeshBlock* m_block1 = mm->createBlock(mbbi);
  
  mm->endCreate(subDomain()->isContinue());
  
  Integer nb_cell = allCells().size();
  Integer nb_env = m_block1->nbEnvironment();
  m_nb_env = nb_env;
  m_nb_vars_to_project = 3 * nb_env + 3 + 1 + 1;
  m_sens_projection = 0;
  // (volumes, ma;sses, energies internes) * nbmatmax
  // + vitesses + energie cinétique + pseudo*
  //
  // on redimensionne les tableaux de la projection en fonction
  // du nombre total d'environnements
  m_is_dir_face.resize(3); // dimension 3
  m_outer_face_normal.resize(6); // dimension 6 faces par mailles
  
  options()->remap()->resizeRemapVariables(m_nb_vars_to_project, m_nb_env);
  
  Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
   ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real3 somme = {0. , 0. , 0.};
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      somme += m_node_coord[inode];
    m_cell_coord[cell] = one_over_nbnode * somme ;
  }
  
  IMeshMaterial* m_mat[nb_env];
  IMeshEnvironment* m_env[nb_env];
  info() << " Initialisation du cas test";
  
  options()->casModel()->initMat(m_dimension);
      
  for (Integer i= 0 ; i <  nb_env; ++i) {
    m_mat[i] = mm->environments()[i]->materials()[0];
    m_env[i] = mm->environments()[0];
  }
  
  
  Int32UniqueArray mat_indexes[nb_env];
  info() << " Nb environments " << nb_env ;
  
  info() << " Trie par environnements  ";
  ENUMERATE_CELL(icell,allCells()){
    Cell cell = *icell;
    if (m_materiau[icell] == 0.) {
      mat_indexes[0].add(icell.itemLocalId());
    } else if (m_materiau[icell] == 1.) {
      mat_indexes[1].add(icell.itemLocalId());
    } else {
      mat_indexes[0].add(icell.itemLocalId());
      mat_indexes[1].add(icell.itemLocalId());
    }
        
  }

  info() << " Ajout de mailles mixtes " ;
 
  
  info() << " fin de la boucle sur les mailles  ";
  Materials::MeshMaterialModifier modifier(mm);
  for (Integer i= 0 ; i <  nb_env; ++i) {
   info() << " mat " << i << " : " << mat_indexes[i].size();
   if (m_mat[i]){
     modifier.addCells(m_mat[i], mat_indexes[i]);
     info() << " mat " << i << " modifié ";
   }
  }
 
  
}

/*---------------------------------------------------------------------------*/
/* A appeler après hydroStartInitEnvAndMat pour préparer                     */
/* traitement des environnements sur accélérateur                            */
/*---------------------------------------------------------------------------*/
void MahycoModule::
_initEnvForAcc() {
  debug() << "_initEnvForAcc";

  m_menv_queue = new MultiAsyncRunQueue(m_runner, mm->environments().size());

  // On récupère sur CPU les adiabatic_cst de chaque environnement
  // On utilise un NumArray pour qu'il soit utilisable aussi sur GPU
  m_adiabatic_cst_env.resize(mm->environments().size());
  Span<Real> out_adiabatic_cst_env(m_adiabatic_cst_env.to1DSpan());
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
    out_adiabatic_cst_env[env->id()] = adiabatic_cst;
  }

  // construit le tableau multi-env m_global_cell_id et le tableau global m_env_id
  _computeMultiEnvGlobalCellId();

  _prepareEnvForAcc();
}

/*---------------------------------------------------------------------------*/
/* Préparer les données multi-envronnement pour l'accélérateur               */
/* A appeler quand la carte des environnements change                        */
/*---------------------------------------------------------------------------*/
void MahycoModule::
_prepareEnvForAcc() {
#ifdef ARCANE_HAS_CUDA
  // "Conseils" utilisation de la mémoire unifiée
  int device = -1;
  cudaGetDevice(&device);

  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    mem_adv_set_read_mostly(env->pureEnvItems().valueIndexes(), device);
  }
#endif
}

/**
 *******************************************************************************
 * \file PrepareFaceGroup()
 * \brief Creation des groupes de faces suivant X, Y et Z
 *******************************************************************************
 */
void MahycoModule::PrepareFaceGroup() {
  Int32UniqueArray face_x0_lid;
  Int32UniqueArray face_y0_lid;
  Int32UniqueArray face_z0_lid;
  Int32UniqueArray face_xmax_lid;
  Int32UniqueArray face_ymax_lid;
  Int32UniqueArray face_zmax_lid;
  Real3 ex = {1. , 0. , 0.};
  Real3 ey = {0. , 1. , 0.};
  Real3 ez = {0. , 0. , 1.};
  Real3 maxCoor= {-1. , -1. , -1.};
  ENUMERATE_NODE(inode, allNodes()){
      maxCoor.x = std::max(maxCoor.x, m_node_coord[inode].x);
      maxCoor.y = std::max(maxCoor.y, m_node_coord[inode].y);
      maxCoor.z = std::max(maxCoor.z, m_node_coord[inode].z);
  }
  ENUMERATE_FACE(iface, allFaces()){
     const Face& face = *iface;
     Integer face_local_id = face.localId();
     for (Integer idir = 0 ; idir <  mesh()->dimension() ; ++idir) {
         m_is_dir_face[face][idir] = false;
     }
     bool flag_x0(true);
     bool flag_y0(true);
     bool flag_z0(true);
     bool flag_xmax(true);
     bool flag_ymax(true);
     bool flag_zmax(true);
     for (NodeEnumerator inode(face.nodes()); inode.index() < face.nbNode(); ++inode) {
         if (m_node_coord[inode].x > options()->threshold)  flag_x0 = false;
         if (m_node_coord[inode].y > options()->threshold)  flag_y0 = false;
         if (m_node_coord[inode].z > options()->threshold)  flag_z0 = false;
         if (math::abs(m_node_coord[inode].x - maxCoor.x) > options()->threshold) flag_xmax = false;
         if (math::abs(m_node_coord[inode].y - maxCoor.y) > options()->threshold) flag_ymax = false;
         if (math::abs(m_node_coord[inode].z - maxCoor.z) > options()->threshold) flag_zmax = false;
     }
     if (flag_x0 == true) face_x0_lid.add(face_local_id);
     if (flag_y0 == true) face_y0_lid.add(face_local_id);
     if (flag_z0 == true) face_z0_lid.add(face_local_id);
     if (flag_xmax == true) face_xmax_lid.add(face_local_id);
     if (flag_ymax == true) face_ymax_lid.add(face_local_id);
     if (flag_zmax == true) face_zmax_lid.add(face_local_id);
   }
   
   mesh()->faceFamily()->createGroup("XMIN", face_x0_lid,true);
   mesh()->faceFamily()->createGroup("YMIN", face_y0_lid,true);
   mesh()->faceFamily()->createGroup("ZMIN", face_z0_lid,true);
   FaceGroup facexmin = mesh()->faceFamily()->findGroup("XMIN");
   FaceGroup faceymin = mesh()->faceFamily()->findGroup("YMIN");
   FaceGroup facezmin = mesh()->faceFamily()->findGroup("ZMIN");
   info() << " taille x 0 " << facexmin.size();
   info() << " taille y 0 " << faceymin.size();
   info() << " taille z 0 " << facezmin.size();
   
   mesh()->faceFamily()->createGroup("XMAX", face_xmax_lid,true);
   mesh()->faceFamily()->createGroup("YMAX", face_ymax_lid,true);
   mesh()->faceFamily()->createGroup("ZMAX", face_zmax_lid,true);
   FaceGroup facexmax = mesh()->faceFamily()->findGroup("XMAX");
   FaceGroup faceymax = mesh()->faceFamily()->findGroup("YMAX");
   FaceGroup facezmax = mesh()->faceFamily()->findGroup("ZMAX");
   info() << " taille x max " << facexmax.size();
   info() << " taille y max " << faceymax.size();
   info() << " taille z max " << facezmax.size();
   
 
   info() << " nombre total de face " << allFaces().size();
   
   info() << " creation des groupes de dimension " << m_dimension;
} 

/*---------------------------------------------------------------------------*/
/* Les listes de faces XMIN, XMAX, YMIN ... doivent être construites au      */
/* préalable par un appel à PrepareFaceGroup()                               */
/*---------------------------------------------------------------------------*/
void MahycoModule::
_initBoundaryConditionsForAcc() {
  debug() << "_initBoundaryConditionsForAcc";
 
  // Remplit la structure contenant les informations sur les conditions aux limites
  // Cela permet de garantir avec les accélérateurs qu'on pourra accéder
  // de manière concurrente aux données.
  {
    m_boundary_conditions.clear();
    for (Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i){
      String NomBC = options()->boundaryCondition[i]->surface;
      FaceGroup face_group = mesh()->faceFamily()->findGroup(NomBC);
      Real value = options()->boundaryCondition[i]->value();
      TypesMahyco::eBoundaryCondition type = options()->boundaryCondition[i]->type();

      BoundaryCondition bcn;
      bcn.nodes = face_group.nodeGroup();
      // attention, cette vue sera à reconstruire si bcn.nodes est modifié
      bcn.boundary_nodes = bcn.nodes.view(); 
      bcn.value = value;
      bcn.type = type;
      m_boundary_conditions.add(bcn);
    }
  }
}

/*---------------------------------------------------------------------------*/
/* Calcul des cell_id globaux : permet d'associer à chaque maille impure (mixte) */
/* l'identifiant de la maille globale                                        */
/*---------------------------------------------------------------------------*/
void MahycoModule::
_computeMultiEnvGlobalCellId() {
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "_computeMultiEnvGlobalCellId";

  // Calcul des cell_id globaux 
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Integer cell_id = cell.localId();
    m_global_cell[cell] = cell_id;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;
        m_global_cell[ev] = cell_id;
      }
      // Maille mixte, contient l'opposé du nombre d'environnements
      m_env_id[icell] = -all_env_cell.nbEnvironment();
    } else {
      // Maille pure, cette boucle est de taille 1
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;
        // Cette affectation n'aura lieu qu'une fois
        m_env_id[icell] = ev.environmentId();
      }
    }
  }

  _checkMultiEnvGlobalCellId();
  PROF_ACC_END;
}

void MahycoModule::
_checkMultiEnvGlobalCellId() {
#ifdef ARCANE_DEBUG
  debug() << "_checkMultiEnvGlobalCellId";

  // Vérification
  ENUMERATE_ENV(ienv, mm) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell(ev.globalCell());
      ARCANE_ASSERT(cell.localId()==m_global_cell[ev], ("lid differents"));
      AllEnvCell all_env_cell(ev.allEnvCell());
      if (all_env_cell.nbEnvironment()==1) {
        ARCANE_ASSERT(m_env_id[cell]==env_id, ("cell pure : environnement id incorrect dans m_env_id[cell]"));
      } else {
        ARCANE_ASSERT(m_env_id[cell]==-all_env_cell.nbEnvironment(), ("cell mixte : m_env_id[cell] différent de -nbEnvironment()"));
      }
    }
  }
#endif
}

