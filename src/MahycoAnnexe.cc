// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"


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
    if (i<3){
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
  // m_nb_vars_to_project = nbmatmax *
  // (0-volumes, 1-masses, 2-energies internes, 
  //  -3,4,5,6- phases, 7-pseudo, 
  // -8,9,10,11,12- deviateurs , -13,14-deformations plastiques) 
  m_nb_vars_to_project = 15 * nb_env ;
  m_sens_projection = 0;
  
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
  if (options()->casModel()->isInternalModel() == true) { 
    info() << " Initialisation du cas test interne au code ";
    options()->casModel()->initMat(m_dimension);
  }
  
  for (Integer i= 0 ; i <  nb_env; ++i) {
    m_mat[i] = mm->environments()[i]->materials()[0];
    m_env[i] = mm->environments()[0];
  }
  
  
  Int32UniqueArray mat_indexes[nb_env];
  info() << " Nb environments " << nb_env ;
  info() << " Nb cells " << allCells().size();
  
  info() << " Trie par environnements  ";
  ENUMERATE_CELL(icell,allCells()){
    Cell cell = *icell;
    if (m_materiau[icell] == 0.) {
      mat_indexes[0].add(icell.itemLocalId());
    } else if (m_materiau[icell] == 1.) {
      mat_indexes[1].add(icell.itemLocalId());
    } else if (m_materiau[icell] == 2.) {
      mat_indexes[2].add(icell.itemLocalId());
    } else if (m_materiau[icell] > 0. && m_materiau[icell] < 1.) {
      mat_indexes[0].add(icell.itemLocalId());
      mat_indexes[1].add(icell.itemLocalId());
    } else if (m_materiau[icell] > 1. && m_materiau[icell] < 2.) {
      mat_indexes[1].add(icell.itemLocalId());
      mat_indexes[2].add(icell.itemLocalId());
    } else {
      mat_indexes[0].add(icell.itemLocalId());
      mat_indexes[1].add(icell.itemLocalId());
      mat_indexes[2].add(icell.itemLocalId());
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
/**
 *******************************************************************************
 * \file PrepareFaceGroup()
 * \brief Creation des groupes de faces suivant X, Y et Z
 *******************************************************************************
 */
void MahycoModule::PrepareFaceGroup() {
    
  if (options()->casModel()->isInternalModel() == true) { 
    Int32UniqueArray face_x0_lid;
    Int32UniqueArray face_y0_lid;
    Int32UniqueArray face_z0_lid;
    Int32UniqueArray face_xmax_lid;
    Int32UniqueArray face_ymax_lid;
    Int32UniqueArray face_zmax_lid;
    Real3 ex = {1. , 0. , 0.};
    Real3 ey = {0. , 1. , 0.};
    Real3 ez = {0. , 0. , 1.};
    Real3 maxCoor= {0. , 0. , 0.};
    Real3 minCoor= {100. , 100. , 100.};
    ENUMERATE_NODE(inode, allNodes()){
        maxCoor.x = std::max(maxCoor.x, m_node_coord[inode].x);
        maxCoor.y = std::max(maxCoor.y, m_node_coord[inode].y);
        maxCoor.z = std::max(maxCoor.z, m_node_coord[inode].z);
        minCoor.x = std::min(minCoor.x, m_node_coord[inode].x);
        minCoor.y = std::min(minCoor.y, m_node_coord[inode].y);
        minCoor.z = std::min(minCoor.z, m_node_coord[inode].z);
    }
    maxCoor.x = parallelMng()->reduce(Parallel::ReduceMax, maxCoor.x);
    maxCoor.y = parallelMng()->reduce(Parallel::ReduceMax, maxCoor.y);
    maxCoor.z = parallelMng()->reduce(Parallel::ReduceMax, maxCoor.z);
    minCoor.x = parallelMng()->reduce(Parallel::ReduceMin, minCoor.x);
    minCoor.y = parallelMng()->reduce(Parallel::ReduceMin, minCoor.y);
    minCoor.z = parallelMng()->reduce(Parallel::ReduceMin, minCoor.z);
    
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
            if (math::abs(m_node_coord[inode].x - minCoor.x) > options()->threshold)  flag_x0 = false;
            if (math::abs(m_node_coord[inode].y - minCoor.y) > options()->threshold)  flag_y0 = false;
            if (math::abs(m_node_coord[inode].z - minCoor.z) > options()->threshold)  flag_z0 = false;
            if (math::abs(m_node_coord[inode].x - maxCoor.x) > options()->threshold) flag_xmax = false;
            if (math::abs(m_node_coord[inode].y - maxCoor.y) > options()->threshold) flag_ymax = false;
            if (math::abs(m_node_coord[inode].z - maxCoor.z) > options()->threshold) flag_zmax = false;
            //if (maxCoor.x > m_node_coord[inode].x) flag_xmax = false;
            //if (maxCoor.y > m_node_coord[inode].y) flag_ymax = false;
            //if (maxCoor.z > m_node_coord[inode].z) flag_zmax = false;
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
    pinfo() << " taille x 0 " << facexmin.size();
    pinfo() << " taille y 0 " << faceymin.size();
    pinfo() << " taille z 0 " << facezmin.size();
    pinfo() << " thresold " << options()->threshold;
    
    mesh()->faceFamily()->createGroup("XMAX", face_xmax_lid,true);
    mesh()->faceFamily()->createGroup("YMAX", face_ymax_lid,true);
    mesh()->faceFamily()->createGroup("ZMAX", face_zmax_lid,true);
    FaceGroup facexmax = mesh()->faceFamily()->findGroup("XMAX");
    FaceGroup faceymax = mesh()->faceFamily()->findGroup("YMAX");
    FaceGroup facezmax = mesh()->faceFamily()->findGroup("ZMAX");
    pinfo() << " taille x max " << facexmax.size();
    pinfo() << " taille y max " << faceymax.size();
    pinfo() << " taille z max " << facezmax.size();
    
    
    info() << " nombre total de face " << allFaces().size();
    
    info() << " creation des groupes de dimension " << m_dimension;
  } 
} 
/**
 *******************************************************************************
 * \file SortieHistory()
 * \brief Creation des groupes de faces suivant X, Y et Z
 *******************************************************************************
 */
void MahycoModule::SortieHistory() {
  if (options()->timeHistory.size() == 0) return;
  
  
  std::ofstream fichier("time-history.txt", std::ofstream::app );
  std::ofstream fichier_noeud("time-history_noeud.txt", std::ofstream::app );
  for (Integer i = 0, nb = options()->timeHistory.size(); i < nb; ++i){
    // Integer period=options()->outputHistoryPeriod();
    Integer period = options()->timeHistory[i]->periode;
    Real frequence = options()->timeHistory[i]->frequence; // not used yet
    if ((m_global_iteration()%period ==0)) { 
      if (fichier.is_open()) { 
        ENUMERATE_ENV(ienv,mm){
        IMeshEnvironment* env = *ienv;
            ENUMERATE_ENVCELL(ienvcell,env)
            {
                EnvCell ev = *ienvcell;   
                Cell cell = ev.globalCell();
                
                if ( m_maille_th() == cell.localId()) {
                    fichier << " Temps " << m_global_time() << " ";
                    fichier << m_maille_th() << " ";
                    // fichier << m_cell_coord[cell] << " ";
                    fichier << m_density[ev] << " ";
                    fichier << m_internal_energy[ev] << " ";
                    fichier << m_pressure[ev] << " ";
                    fichier << m_temperature[ev] << " ";
                    fichier << m_sound_speed[ev] << " ";
                    fichier << m_frac_phase1[ev] << " ";
                    fichier << m_frac_phase2[ev] << " ";
                    fichier << m_frac_phase3[ev] << " ";
                    fichier << m_strain_norm[ev] << " ";
                    fichier << m_strain_tensor_xx[ev] << " ";
                    fichier << m_strain_tensor_yy[ev] << " ";
                    fichier << m_strain_tensor_xy[ev] << " ";
                    fichier << m_strain_tensor_yz[ev] << " ";
                    fichier << m_strain_tensor_xz[ev] << " ";
                    fichier << std::endl;
                }
             }
          }
      }
      if (fichier_noeud.is_open()) { 
        ENUMERATE_NODE(inode, allNodes()){
            Node node = *inode;
            if ( m_noeud_th() == node.localId()) {
                fichier_noeud << " Temps " << m_global_time() << " ";
                fichier_noeud << m_noeud_th() << " ";
                fichier_noeud << m_node_coord[node].x << " ";
                fichier_noeud << m_node_coord[node].y << " ";
                fichier_noeud << m_node_coord[node].z << " ";
                fichier_noeud << m_velocity[node].x << " ";
                fichier_noeud << m_velocity[node].y << " ";
                fichier_noeud << m_velocity[node].z << " ";
                fichier_noeud << std::endl;
            }
        }
      }
    }
  }
}
/**
 *******************************************************************************
 * \file lireFichierPourCDL()
 * \brief 
 *******************************************************************************
 */
void MahycoModule::lireFichierCDL(const std::string& nomFichier)
{
    std::ifstream fichier(nomFichier); // Ouverture du fichier en lecture
    pinfo() << "Ouverture de " << nomFichier;
    
    if (fichier)
    {
        
        int compteur = 0;
        Real x(-1),y(-1);

        while (fichier >> x >> y)
        {
            compteur++;
            pinfo() << x << " " << y;
            m_table_temps.push_back(x);
            m_table_pression.push_back(y);
        }
        m_taille_table = compteur;
        
        pinfo() << " Taille de la table de marche des CDL de pression" << m_taille_table;
        fichier.close(); // Fermeture du fichier
    }
    else
    {
        std::cout << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
}
/**
 *******************************************************************************
 * \file InitTH()
 * \brief 
 *******************************************************************************
 */
void MahycoModule::initTH()
{
  // pour l'instant : options()->timeHistory.size() est de taille 1
  for (Integer i = 0, nb = options()->timeHistory.size(); i < nb; ++i){
    Real3 borneSup = options()->timeHistory[i]->borneSup; 
    Real3 borneInf = options()->timeHistory[i]->borneInf; 
    // Sortie Time History
    m_maille_th = -1;
    bool trouve = false;
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      if ((borneInf.x < m_cell_coord[cell].x) && (m_cell_coord[cell].x < borneSup.x) 
        && (borneInf.y < m_cell_coord[cell].y) && (m_cell_coord[cell].y < borneSup.y)
        && (borneInf.z < m_cell_coord[cell].z) && (m_cell_coord[cell].z < borneSup.z)
        && cell.isOwn() ) {
        m_maille_th = cell.localId();
        trouve = true;
      }
    }
    if (!trouve) {
        std::cout << "Pas de maille trouvé pour le time history dans ce sous-domaine ( m_maille_th = " << m_maille_th() << " )" << std::endl;
    } else {
        std::cout << " Maille " << m_maille_th() << " trouvé dans le sous-domaine " << my_rank << std::endl;
    }
    Real3 borneSupNoeud = options()->timeHistory[i]->borneSupNoeud; 
    Real3 borneInfNoeud = options()->timeHistory[i]->borneInfNoeud; 
    // Sortie Time History pour les noeuds
    m_noeud_th = -1;
    ENUMERATE_NODE(inode, allNodes()){
        Node node = *inode;
         if ((borneInfNoeud.x < m_node_coord[node].x) && (m_node_coord[node].x < borneSupNoeud.x) 
        && (borneInfNoeud.y < m_node_coord[node].y) && (m_node_coord[node].y < borneSupNoeud.y)
        && (borneInfNoeud.z < m_node_coord[node].z) && (m_node_coord[node].z < borneSupNoeud.z)
        && node.isOwn() ) {
        m_noeud_th = node.localId();
        trouve = true;
      }
    }
    if (!trouve) {
        std::cout << "Pas de noeud trouvé pour le time history dans ce sous-domaine ( m_noeud_th = " << m_noeud_th() << " )" << std::endl;
    } else {
        std::cout << " Noeud " << m_noeud_th() << " trouvé dans le sous-domaine " << my_rank << std::endl;
    }
    
    
  }
}
