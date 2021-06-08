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
  IMeshMaterial* m_mat1;
  IMeshMaterial* m_mat2;
  IMeshEnvironment* m_env1;
  IMeshEnvironment* m_env2;

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
  
   ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real3 somme = {0. , 0. , 0.};
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      somme += m_node_coord[inode];
    m_cell_coord[cell] = 0.125 * somme ;
  }
  
  info() << " Initialisation du cas test";
  
  options()->casModel()->initMat();
      
  m_mat1 = mm->environments()[0]->materials()[0];
  m_mat2 = mm->environments()[1]->materials()[0];
  m_env1 = mm->environments()[0];
  m_env2 = mm->environments()[1];
  
  
  Int32UniqueArray mat1_indexes;
  Int32UniqueArray mat2_indexes;
  info() << " Nb environments " << nb_env ;
  
  info() << " Trie par environnements  ";
  ENUMERATE_CELL(icell,allCells()){
    Cell cell = *icell;
    Int64 cell_index = cell.uniqueId();
    if (m_materiau[icell] == 0.) {
      mat1_indexes.add(icell.itemLocalId());
    } else if (m_materiau[icell] == 1.) {
      mat2_indexes.add(icell.itemLocalId());
    } else {
      mat1_indexes.add(icell.itemLocalId());
      mat2_indexes.add(icell.itemLocalId());
    }
        
  }

  info() << " Ajout de mailles mixtes " ;
 
  
  info() << " fin de la boucle sur les mailles  ";
  info() << " mat 1 " << mat1_indexes.size();
  info() << " mat 2 " << mat2_indexes.size();
  if (m_mat1){
    Materials::MeshMaterialModifier modifier(mm);
    modifier.addCells(m_mat1, mat1_indexes);
    info() << " mat1 modifie ";
  }    
  if (m_mat2){
    Materials::MeshMaterialModifier modifier(mm);
    modifier.addCells(m_mat2, mat2_indexes);
    info() << " mat2 modifie ";
  }
  
}
/**
 *******************************************************************************
 * \file PrepareFaceGroup()
 * \brief Creation des groupes de faces suivant X, Y et Z
 *******************************************************************************
 */
void MahycoModule::PrepareFaceGroup() {
  Int32UniqueArray face_x_lid;
  Int32UniqueArray face_y_lid;
  Int32UniqueArray face_z_lid;
  Int32UniqueArray face_x0_lid;
  Int32UniqueArray face_y0_lid;
  Int32UniqueArray face_z0_lid;
  Int32UniqueArray face_xmax_lid;
  Int32UniqueArray face_ymax_lid;
  Int32UniqueArray face_zmax_lid;
  Int32UniqueArray face_interne_x_lid;
  Int32UniqueArray face_interne_y_lid;
  Int32UniqueArray face_interne_z_lid;
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
     for (Integer idir = 0 ; idir <  3 ; ++idir) {
         m_is_dir_face[face][idir] = false;
     }
     // on prend un vecteur quelconque de la 
     Real3 face_vec1 = m_node_coord[face.node(2)] - m_node_coord[face.node(0)]; 
     Real3 face_vec2 = m_node_coord[face.node(3)] - m_node_coord[face.node(1)]; 
     Real productx = math::dot(face_vec2, ex);
     Real producty = math::dot(face_vec2, ey);
     Real productz = math::dot(face_vec1, ez);
     if (m_node_coord[face.node(0)].x < options()->threshold 
         && m_node_coord[face.node(1)].x < options()->threshold 
         && m_node_coord[face.node(2)].x < options()->threshold 
         && m_node_coord[face.node(3)].x < options()->threshold) {
         face_x0_lid.add(face_local_id);
     } else if (m_node_coord[face.node(0)].y < options()->threshold 
         && m_node_coord[face.node(1)].y < options()->threshold 
         && m_node_coord[face.node(2)].y < options()->threshold 
         && m_node_coord[face.node(3)].y < options()->threshold) {
         face_y0_lid.add(face_local_id);
     } else if (m_node_coord[face.node(0)].z < options()->threshold 
         && m_node_coord[face.node(1)].z < options()->threshold 
         && m_node_coord[face.node(2)].z < options()->threshold 
         && m_node_coord[face.node(3)].z < options()->threshold) {
         face_z0_lid.add(face_local_id);
     }
     if (math::abs(m_node_coord[face.node(0)].x - maxCoor.x)< options()->threshold 
         && math::abs(m_node_coord[face.node(1)].x - maxCoor.x)< options()->threshold 
         && math::abs(m_node_coord[face.node(2)].x - maxCoor.x)< options()->threshold 
         && math::abs(m_node_coord[face.node(3)].x - maxCoor.x)< options()->threshold) {
         face_xmax_lid.add(face_local_id);
     } else if (math::abs(m_node_coord[face.node(0)].y - maxCoor.y)< options()->threshold 
         && math::abs(m_node_coord[face.node(1)].y - maxCoor.y)< options()->threshold 
         && math::abs(m_node_coord[face.node(2)].y - maxCoor.y)< options()->threshold 
         && math::abs(m_node_coord[face.node(3)].y - maxCoor.y)< options()->threshold) {
         face_ymax_lid.add(face_local_id);
     } else if (math::abs(m_node_coord[face.node(0)].z - maxCoor.z) < options()->threshold 
         && math::abs(m_node_coord[face.node(1)].z - maxCoor.z)< options()->threshold 
         && math::abs(m_node_coord[face.node(2)].z - maxCoor.z)< options()->threshold 
         && math::abs(m_node_coord[face.node(3)].z - maxCoor.z)< options()->threshold) {
         face_zmax_lid.add(face_local_id);
     }
     if (math::abs(productx) < options()->threshold) {
         face_x_lid.add(face_local_id);
         m_is_dir_face[face][0] = true;
         if (face.nbCell() == 2) face_interne_x_lid.add(face_local_id);
     } else if (math::abs(producty) < options()->threshold) {
         face_y_lid.add(face_local_id);
         m_is_dir_face[face][1] = true;
         if (face.nbCell() == 2) face_interne_y_lid.add(face_local_id);
     } else if (math::abs(productz) < options()->threshold) {
         face_z_lid.add(face_local_id);
         m_is_dir_face[face][2] = true;
         if (face.nbCell() == 2) face_interne_z_lid.add(face_local_id);
     }
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
   
   mesh()->faceFamily()->createGroup("FACE_X", face_x_lid,true);
   mesh()->faceFamily()->createGroup("FACE_X_INTERNE", face_interne_x_lid,true);
   FaceGroup inner_x_faces = mesh()->faceFamily()->findGroup("FACE_X");
   FaceGroup inner_interne_x_faces = mesh()->faceFamily()->findGroup("FACE_X_INTERNE");
   info() << " taille x interne " << inner_interne_x_faces.size();
   info() << " taille x " << inner_x_faces.size();
   
   mesh()->faceFamily()->createGroup("FACE_Y", face_y_lid,true);
   mesh()->faceFamily()->createGroup("FACE_Y_INTERNE", face_interne_y_lid,true);
   FaceGroup inner_y_faces = mesh()->faceFamily()->findGroup("FACE_Y");
   FaceGroup inner_interne_y_faces = mesh()->faceFamily()->findGroup("FACE_Y_INTERNE");
   info() << " taille y interne " << inner_interne_y_faces.size();
   info() << " taille y " << inner_y_faces.size();
   
   mesh()->faceFamily()->createGroup("FACE_Z", face_z_lid,true);
   mesh()->faceFamily()->createGroup("FACE_Z_INTERNE", face_interne_z_lid,true);
   FaceGroup inner_z_faces = mesh()->faceFamily()->findGroup("FACE_Z");
   FaceGroup inner_interne_z_faces = mesh()->faceFamily()->findGroup("FACE_Z_INTERNE");
   info() << " taille z interne " << inner_interne_z_faces.size();
   info() << " taille z " << inner_z_faces.size();
   
   info() << " taille totale confirmee " << inner_x_faces.size() + inner_y_faces.size() + inner_z_faces.size();
   info() << " taille totale interne confirmee " << inner_interne_x_faces.size() + inner_interne_y_faces.size() + inner_interne_z_faces.size();
   info() << " nombre total de face " << allFaces().size();
   
   info() << " creation des groupes OK ";
} 
