﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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
    // Le bloc ne contient que 5 milieux au max
    if (i<5){
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
  // 15-energies internes old
  // 16-temperature
  // 17-temperature-old
  m_nb_vars_to_project = 18 * nb_env ;
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
    } else if (m_materiau[icell] > 0. && m_materiau[icell] < 1.) {
      mat_indexes[0].add(icell.itemLocalId());
      mat_indexes[1].add(icell.itemLocalId());
    } else if (m_materiau[icell] == 1.) {
      mat_indexes[1].add(icell.itemLocalId());
    } else if (m_materiau[icell] > 1. && m_materiau[icell] < 2.) {
      mat_indexes[1].add(icell.itemLocalId());
      mat_indexes[2].add(icell.itemLocalId());
    } else if (m_materiau[icell] == 2.) {
      mat_indexes[2].add(icell.itemLocalId());
    } else if (m_materiau[icell] > 2. && m_materiau[icell] < 3.) {
      mat_indexes[2].add(icell.itemLocalId());
      mat_indexes[3].add(icell.itemLocalId());
    } else if (m_materiau[icell] == 3.) {
      mat_indexes[3].add(icell.itemLocalId());
    } else {
      mat_indexes[0].add(icell.itemLocalId());
      mat_indexes[1].add(icell.itemLocalId());
      mat_indexes[2].add(icell.itemLocalId());
      mat_indexes[3].add(icell.itemLocalId());
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
    // pinfo() << " taille x 0 " << facexmin.size();
    // pinfo() << " taille y 0 " << faceymin.size();
    // pinfo() << " taille z 0 " << facezmin.size();
    // pinfo() << " thresold " << options()->threshold;
    
    mesh()->faceFamily()->createGroup("XMAX", face_xmax_lid,true);
    mesh()->faceFamily()->createGroup("YMAX", face_ymax_lid,true);
    mesh()->faceFamily()->createGroup("ZMAX", face_zmax_lid,true);
    FaceGroup facexmax = mesh()->faceFamily()->findGroup("XMAX");
    FaceGroup faceymax = mesh()->faceFamily()->findGroup("YMAX");
    FaceGroup facezmax = mesh()->faceFamily()->findGroup("ZMAX");
    // pinfo() << " taille x max " << facexmax.size();
    // pinfo() << " taille y max " << faceymax.size();
    // pinfo() << " taille z max " << facezmax.size();
    
    
    // info() << " nombre total de face " << allFaces().size();
    
    // info() << " creation des groupes de dimension " << m_dimension;
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
