// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* CartesianMesh.cc                                            (C) 2000-2014 */
/*                                                                           */
/* Maillage cartésien.                                                       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/ArcanePrecomp.h"

#include "arcane/utils/TraceAccessor.h"
#include "arcane/utils/NotImplementedException.h"
#include "arcane/utils/AutoDestroyUserData.h"
#include "arcane/utils/IUserDataList.h"

#include "arcane/IMesh.h"
#include "arcane/ItemPrinter.h"
#include "arcane/IItemFamily.h"
#include "arcane/IParallelMng.h"
#include "arcane/VariableTypes.h"
#include "arcane/Properties.h"

#include "cartesian/ICartesianMesh.h"
#include "cartesian/CellDirectionMng.h"
#include "cartesian/FaceDirectionMng.h"
#include "cartesian/NodeDirectionMng.h"
#include "cartesian/CartesianConnectivity.h"
#include "cartesian/FactCartDirectionMng.h"
#include "cartesian/CartesianGridT.h"
#include "cartesian/NumberingConverterT.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Infos spécifiques à un maillage cartésien.
 */
class CartesianMesh
: public TraceAccessor
, public ICartesianMesh
{
 public:
  //! Type d'une grille cartésienne en numérotation locale
  using CartesianGrid = CartesianGridT<LocalIdType>;
  //! Type pour convertir numérotation Face => Cell
  using NumbConvFaceCell = NumberingConverterT<Face,Cell>;

 public:
  CartesianMesh(IMesh* mesh)
  : TraceAccessor(mesh->traceMng()), m_mesh(mesh),
  m_fact_cart_dm (m_mesh),
  m_cart_grid (*m_fact_cart_dm.cartesianGrid()),
  m_connectivity(mesh->itemsInternal(IK_Cell).data(), mesh->itemsInternal(IK_Node).data(), m_cart_grid)
  {
    Integer nb_dir = mesh->dimension();
    for( Integer i=0; i<nb_dir; ++i ){
      m_local_face_direction[i] = -1;
      m_cell_directions[i].init(this,(eMeshDirection)i, m_fact_cart_dm.cartesianGrid()->cartNumCellPtr());
      m_numb_conv_face_cell[i] = new NumbConvFaceCell(i, m_cart_grid);
      m_face_directions[i].init(this,(eMeshDirection)i, m_fact_cart_dm.cartesianGrid(), m_numb_conv_face_cell[i]);
      m_node_directions[i].init(this,(eMeshDirection)i, m_fact_cart_dm.cartesianGrid()->cartNumNodePtr());
    }
  }
  ~CartesianMesh()
  {
    for( Integer i=0; i<3; ++i ){
      m_cell_directions[i].destroy();
      m_face_directions[i].destroy();
      delete m_numb_conv_face_cell[i];
      m_node_directions[i].destroy();
    }
  }
  void build() override {}

  //! Maillage associé à ce maillage cartésien
  IMesh* mesh() const override { return m_mesh; }

  //! Gestionnaire de trace associé.
  ITraceMng* traceMng() const override { return TraceAccessor::traceMng(); }

  CellDirectionMng cellDirection(eMeshDirection dir) override
  {
    return m_cell_directions[dir];
  }

  CellDirectionMng cellDirection(Integer idir) override
  {
    return m_cell_directions[idir];
  }

  FaceDirectionMng faceDirection(eMeshDirection dir) override
  {
    return m_face_directions[dir];
  }

  FaceDirectionMng faceDirection(Integer idir) override
  {
    return m_face_directions[idir];
  }

  NodeDirectionMng nodeDirection(eMeshDirection dir) override
  {
    return m_node_directions[dir];
  }

  NodeDirectionMng nodeDirection(Integer idir) override
  {
    return m_node_directions[idir];
  }

  void computeDirections() override;

  CartesianConnectivity connectivity() override
  {
    return m_connectivity;
  }

  void _computeMeshDirection(eMeshDirection dir,VariableCellReal3& cells_center,
                             VariableFaceReal3& faces_center);
  void _computeNodeCellInformations2D(Cell cell0,Real3 cell0_coord,VariableNodeReal3& nodes_coord);
  void _computeNodeCellInformations3D(Cell cell0,Real3 cell0_coord,VariableNodeReal3& nodes_coord);

 private:

  //! Indice dans la numérotation locale de la maille, de la face dans
  // la direction X, Y ou Z
  Int32 m_local_face_direction[3];
  IMesh* m_mesh;
  CellDirectionMng m_cell_directions[3];
  FaceDirectionMng m_face_directions[3];
  NodeDirectionMng m_node_directions[3];
  NumbConvFaceCell* m_numb_conv_face_cell[3]={nullptr,nullptr,nullptr};
  FactCartDirectionMng m_fact_cart_dm;  //! Fabrique de grille cartésienne
  const CartesianGrid& m_cart_grid;  //! Référence sur la grille cartésienne
  CartesianConnectivity m_connectivity;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

extern "C++" ICartesianMesh*
arcaneCreateCartesianMesh(IMesh* mesh)
{
  CartesianMesh* cm = new CartesianMesh(mesh);
  cm->build();
  return cm;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CartesianMesh::
computeDirections()
{
  VariableCellReal3 cells_center(VariableBuildInfo(m_mesh,"TemporaryCartesianMeshCellCenter"));
  VariableFaceReal3 faces_center(VariableBuildInfo(m_mesh,"TemporaryCartesianMeshFaceCenter"));

  // Calcule les coordonnées du centre des mailles.
  VariableNodeReal3& nodes_coord = m_mesh->nodesCoordinates();
  ENUMERATE_CELL(icell,m_mesh->allCells()){
    Cell cell = *icell;
    Real3 center;
    for( NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode )
      center += nodes_coord[inode];
    center /= cell.nbNode();
    cells_center[icell] = center;
  }
  ENUMERATE_FACE(iface,m_mesh->allFaces()){
    Face face = *iface;
    Real3 center;
    for( NodeEnumerator inode(face.nodes()); inode.hasNext(); ++inode )
      center += nodes_coord[inode];
    center /= face.nbNode();
    faces_center[iface] = center;
  }

  IItemFamily* cell_family = m_mesh->cellFamily();
  Int32 next_face_x = -1;
  Int32 next_face_y = -1;
  Int32 next_face_z = -1;

  CellVectorView cell_view = cell_family->allItems().view();
  Cell cell0 = cell_view[0];
  Integer nb_face = cell0.nbFace();
  Integer nb_node = cell0.nbNode();
  Real3 cell_center = cells_center[cell0];

  info() << "Cartesian mesh compute directions";

  for( Integer i=0; i<nb_node; ++i ){
    Node node = cell0.node(i);
    info() << "Node I=" << i << " node=" << ItemPrinter(node) << " pos=" << nodes_coord[node];
  }

  // On suppose que toutes les mailles ont le même sens de numérotation dans le maillage.
  // Par exemple, pour toutes les mailles, la face d'indice 0 est celle du haut, celle
  // d'indice 1 celle de droite.
  for( Integer i=0; i<nb_face; ++i ){
    Face f = cell0.face(i);
    if (f.isSubDomainBoundary())
      continue;
    Cell next_cell = (f.backCell()==cell0) ? f.frontCell() : f.backCell();
    Real3 next_center = cells_center[next_cell];
    info(4) << "NEXT_CELL=" << ItemPrinter(next_cell) << " center=" << next_center \
            << " back=" << f.backCell().uniqueId()
            << " front=" << f.frontCell().uniqueId();

    Real diff_x = math::abs(next_center.x - cell_center.x);
    Real diff_y = math::abs(next_center.y - cell_center.y);
    Real diff_z = math::abs(next_center.z - cell_center.z);
    info(4) << "NEXT_CELL=" << ItemPrinter(next_cell) << " diff=" << Real3(diff_x,diff_y,diff_z);
    //TODO: Verifier qu'il s'agit bien de la maille apres et pas avant.
    // (tenir compte du signe de diff)
    if (diff_x>diff_y && diff_x>diff_z){
      // INC X
      next_face_x = i;
      info() << "Advance in direction X -> " << next_face_x;
    }
    else if (diff_y>diff_x && diff_y>diff_z){
      // INC Y
      next_face_y = i;
      info() << "Advance in direction Y -> " << next_face_y;
    }
    else if (diff_z>diff_x && diff_z>diff_y){
      // INC Z
      next_face_z = i;
      info() << "Advance in direction Z -> " << next_face_z;
    }
    else
      throw FatalErrorException(A_FUNCINFO,"Bad value for next cell");
  }

  bool is_3d = m_mesh->dimension()==3;
  if (is_3d)
    _computeNodeCellInformations3D(cell0,cells_center[cell0],nodes_coord);
  else
    _computeNodeCellInformations2D(cell0,cells_center[cell0],nodes_coord);

  info() << "Informations from IMesh properties:";
  Properties* mesh_properties = mesh()->properties();
  Int64 global_nb_cell_x = mesh_properties->getInt64WithDefault("GlobalNbCellX",-1);
  Int64 global_nb_cell_y = mesh_properties->getInt64WithDefault("GlobalNbCellY",-1);
  Int64 global_nb_cell_z = mesh_properties->getInt64WithDefault("GlobalNbCellZ",-1);
  info() << "GlobalNbCell: "
         << " X=" << global_nb_cell_x
         << " Y=" << global_nb_cell_y
         << " Z=" << global_nb_cell_z;

  Int32 own_nb_cell_x = mesh_properties->getInt32WithDefault("OwnNbCellX",-1);
  Int32 own_nb_cell_y = mesh_properties->getInt32WithDefault("OwnNbCellY",-1);
  Int32 own_nb_cell_z = mesh_properties->getInt32WithDefault("OwnNbCellZ",-1);
  info() << "OwnNbCell: "
         << " X=" << own_nb_cell_x
         << " Y=" << own_nb_cell_y
         << " Z=" << own_nb_cell_z;

  Int32 sub_domain_offset_x = mesh_properties->getInt32WithDefault("SubDomainOffsetX",-1);
  Int32 sub_domain_offset_y = mesh_properties->getInt32WithDefault("SubDomainOffsetY",-1);
  Int32 sub_domain_offset_z = mesh_properties->getInt32WithDefault("SubDomainOffsetZ",-1);
  info() << "SubDomainOffset: "
         << " X=" << sub_domain_offset_x
         << " Y=" << sub_domain_offset_y
         << " Z=" << sub_domain_offset_z;

  Int64 own_cell_offset_x = mesh_properties->getInt64WithDefault("OwnCellOffsetX",-1);
  Int64 own_cell_offset_y = mesh_properties->getInt64WithDefault("OwnCellOffsetY",-1);
  Int64 own_cell_offset_z = mesh_properties->getInt64WithDefault("OwnCellOffsetZ",-1);
  info() << "OwnCellOffset: "
         << " X=" << own_cell_offset_x
         << " Y=" << own_cell_offset_y
         << " Z=" << own_cell_offset_z;

  if (next_face_x!=(-1)){
    m_local_face_direction[MD_DirX] = next_face_x;
    _computeMeshDirection(MD_DirX,cells_center,faces_center);
    CellDirectionMng& cdm = m_cell_directions[MD_DirX];
    cdm.m_global_nb_cell = global_nb_cell_x;
    cdm.m_own_nb_cell = own_nb_cell_x;
    cdm.m_sub_domain_offset = sub_domain_offset_x;
    cdm.m_own_cell_offset = own_cell_offset_x;
  }
  if (next_face_y!=(-1)){
    m_local_face_direction[MD_DirY] = next_face_y;
    _computeMeshDirection(MD_DirY,cells_center,faces_center);
    CellDirectionMng& cdm = m_cell_directions[MD_DirY];
    cdm.m_global_nb_cell = global_nb_cell_y;
    cdm.m_own_nb_cell = own_nb_cell_y;
    cdm.m_sub_domain_offset = sub_domain_offset_y;
    cdm.m_own_cell_offset = own_cell_offset_y;
  }
  if (next_face_z!=(-1)){
    m_local_face_direction[MD_DirZ] = next_face_z;
    _computeMeshDirection(MD_DirZ,cells_center,faces_center);
    CellDirectionMng& cdm = m_cell_directions[MD_DirZ];
    cdm.m_global_nb_cell = global_nb_cell_z;
    cdm.m_own_nb_cell = own_nb_cell_z;
    cdm.m_sub_domain_offset = sub_domain_offset_z;
    cdm.m_own_cell_offset = own_cell_offset_z;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
 * \brief Calcule les infos sur les noeuds avant/après et gauche/droite d'une maille
 * pour chaque direction.
 */
void CartesianMesh::
_computeNodeCellInformations2D(Cell cell0,Real3 cell0_coord,VariableNodeReal3& nodes_coord)
{
  Int32 nodes_indirection_i[CellDirectionMng::MAX_NB_NODE];
  Int32ArrayView nodes_indirection(CellDirectionMng::MAX_NB_NODE,nodes_indirection_i);
  Integer nb_node = cell0.nbNode();
  Real3 cell_coord = cell0_coord;

  // Direction X
  nodes_indirection.fill(-1);
  for( Integer i=0; i<nb_node; ++i ){
    Node node = cell0.node(i);
    Real3 node_coord = nodes_coord[node];
    if (node_coord.x>cell_coord.x){
      if (node_coord.y>cell_coord.y)
        nodes_indirection[CNP_NextLeft] = i;
      else
        nodes_indirection[CNP_NextRight] = i;
    }
    else{
      if (node_coord.y>cell_coord.y)
        nodes_indirection[CNP_PreviousLeft] = i;
      else
        nodes_indirection[CNP_PreviousRight] = i;
    }
  }
  m_cell_directions[MD_DirX].setNodesIndirection(nodes_indirection);

  // Direction Y
  nodes_indirection.fill(-1);
  for( Integer i=0; i<nb_node; ++i ){
    Node node = cell0.node(i);
    Real3 node_coord = nodes_coord[node];
    if (node_coord.y>cell_coord.y){
      if (node_coord.x>cell_coord.x)
        nodes_indirection[CNP_NextRight] = i;
      else
        nodes_indirection[CNP_NextLeft] = i;
    }
    else{
      if (node_coord.x>cell_coord.x)
        nodes_indirection[CNP_PreviousRight] = i;
      else
        nodes_indirection[CNP_PreviousLeft] = i;
    }
  }
  m_cell_directions[MD_DirY].setNodesIndirection(nodes_indirection);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
 * \brief Calcule les infos sur les noeuds avant/après et gauche/droite d'une maille
 * pour chaque direction.
 */
void CartesianMesh::
_computeNodeCellInformations3D(Cell cell0,Real3 cell0_coord,VariableNodeReal3& nodes_coord)
{
  Int32 nodes_indirection_i[CellDirectionMng::MAX_NB_NODE];
  Int32ArrayView nodes_indirection(CellDirectionMng::MAX_NB_NODE,nodes_indirection_i);
  Integer nb_node = cell0.nbNode();
  Real3 cell_coord = cell0_coord;
  bool is_3d = m_mesh->dimension()==3;

  // Direction X
  nodes_indirection.fill(-1);
  for( Integer i=0; i<nb_node; ++i ){
    Node node = cell0.node(i);
    Real3 node_coord = nodes_coord[node];
    if (node_coord.z>cell_coord.z){
      if (node_coord.x>cell_coord.x){
        if (node_coord.y>cell_coord.y)
          nodes_indirection[CNP_TopNextLeft] = i;
        else
          nodes_indirection[CNP_TopNextRight] = i;
      }
      else{
        if (node_coord.y>cell_coord.y)
          nodes_indirection[CNP_TopPreviousLeft] = i;
        else
          nodes_indirection[CNP_TopPreviousRight] = i;
      }
    }
    else{
      if (node_coord.x>cell_coord.x){
        if (node_coord.y>cell_coord.y)
          nodes_indirection[CNP_NextLeft] = i;
        else
          nodes_indirection[CNP_NextRight] = i;
      }
      else{
        if (node_coord.y>cell_coord.y)
          nodes_indirection[CNP_PreviousLeft] = i;
        else
          nodes_indirection[CNP_PreviousRight] = i;
      }
    }
  }
  m_cell_directions[MD_DirX].setNodesIndirection(nodes_indirection);

  // Direction Y
  nodes_indirection.fill(-1);
  for( Integer i=0; i<nb_node; ++i ){
    Node node = cell0.node(i);
    Real3 node_coord = nodes_coord[node];
    if (node_coord.z<cell_coord.z){
      if (node_coord.y>cell_coord.y){
        if (node_coord.x>cell_coord.x)
          nodes_indirection[CNP_TopNextRight] = i;
        else
          nodes_indirection[CNP_TopNextLeft] = i;
      }
      else{
        if (node_coord.x>cell_coord.x)
          nodes_indirection[CNP_TopPreviousRight] = i;
        else
          nodes_indirection[CNP_TopPreviousLeft] = i;
      }
    }
    else{
      if (node_coord.y>cell_coord.y){
        if (node_coord.x>cell_coord.x)
          nodes_indirection[CNP_NextRight] = i;
        else
          nodes_indirection[CNP_NextLeft] = i;
      }
      else{
        if (node_coord.x>cell_coord.x)
          nodes_indirection[CNP_PreviousRight] = i;
        else
          nodes_indirection[CNP_PreviousLeft] = i;
      }
    }
  }
  m_cell_directions[MD_DirY].setNodesIndirection(nodes_indirection);

  if (is_3d){
    nodes_indirection.fill(-1);
    for( Integer i=0; i<nb_node; ++i ){
      Node node = cell0.node(i);
      Real3 node_coord = nodes_coord[node];
      if (node_coord.y>cell_coord.y){
        if (node_coord.z>cell_coord.z){
          if (node_coord.x>cell_coord.x)
            nodes_indirection[CNP_TopNextRight] = i;
          else
            nodes_indirection[CNP_TopNextLeft] = i;
        }
        else{
          if (node_coord.x>cell_coord.x)
            nodes_indirection[CNP_TopPreviousRight] = i;
          else
            nodes_indirection[CNP_TopPreviousLeft] = i;
        }
      }
      else{
        if (node_coord.z>cell_coord.z){
          if (node_coord.x>cell_coord.x)
            nodes_indirection[CNP_NextRight] = i;
          else
            nodes_indirection[CNP_NextLeft] = i;
        }
        else{
          if (node_coord.x>cell_coord.x)
            nodes_indirection[CNP_PreviousRight] = i;
          else
            nodes_indirection[CNP_PreviousLeft] = i;
        }
      }
    }
    m_cell_directions[MD_DirZ].setNodesIndirection(nodes_indirection);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CartesianMesh::
_computeMeshDirection(eMeshDirection dir,VariableCellReal3& /*cells_center*/,
                      VariableFaceReal3& /*faces_center*/)
{
  IItemFamily* cell_family = m_mesh->cellFamily();
  IItemFamily* face_family = m_mesh->faceFamily();
  IItemFamily* node_family = m_mesh->nodeFamily();

  CellDirectionMng& cell_dm = m_cell_directions[dir];
  FaceDirectionMng& face_dm = m_face_directions[dir];
  NodeDirectionMng& node_dm = m_node_directions[dir];

  //TODO: attention a remettre a jour apres changement de maillage.
  //cell_dm.m_items = cell_family->itemsInternal();
  info() << "COMPUTE DIRECTION dir=" << dir;

  Int32 prev_local_face = -1;
  Int32 next_local_face = m_local_face_direction[dir];
  Integer mesh_dim = m_mesh->dimension();
  // Calcul le numero local de face oppose à la face suivante.
  if (mesh_dim==2)
    prev_local_face = (next_local_face + 2) % 4;
  else if (mesh_dim==3)
    prev_local_face = (next_local_face + 3) % 6;

  cell_dm.setLocalFaceIndex(next_local_face,prev_local_face);
  cell_dm.computeInnerAndOuterItems(cell_family->allItems());
  
  bool use_cart_impl=true;
  if (use_cart_impl) {
    face_dm.computeInnerAndOuterItemsFromFamily(face_family);
  } else {
    CellGroup all_cells = cell_family->allItems();
    Int32UniqueArray faces_lid;
    faces_lid.reserve(all_cells.size());

    // Calcule la liste des faces dans une direction donnée.
    ENUMERATE_CELL(icell,all_cells){
      Cell cell = *icell;
      Face next_face = cell.face(next_local_face);
      Face prev_face = cell.face(prev_local_face);

      //! Ajoute la face d'avant à la liste des faces de cette direction
      faces_lid.add(prev_face.localId());
      // Ajoute la face d'après mais uniquement si elle est au bord (car sinon elle a déjà été ajoutée
      // avec l'instruction précédente)
      if (next_face.nbCell()==1)
        faces_lid.add(next_face.localId());


    }
    FaceGroup all_faces_dir = face_family->createGroup(String("AllFacesDirection")+(int)dir,faces_lid,true);

    face_dm.computeInnerAndOuterItems(all_faces_dir);
  }

  node_dm.computeInnerAndOuterItems(node_family->allItems());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ICartesianMesh* ICartesianMesh::
getReference(IMesh* mesh,bool create)
{
  //TODO: faire lock pour multi-thread
  const char* name = "CartesianMesh";
  IUserDataList* udlist = mesh->userDataList();

  IUserData* ud = udlist->data(name,true);
  if (!ud){
    if (!create)
      return 0;
    ICartesianMesh* cm = arcaneCreateCartesianMesh(mesh);
    udlist->setData(name,new AutoDestroyUserData<ICartesianMesh>(cm));
    return cm;
  }
  AutoDestroyUserData<ICartesianMesh>* adud = dynamic_cast<AutoDestroyUserData<ICartesianMesh>*>(ud);
  if (!adud)
    throw FatalErrorException(A_FUNCINFO,"Can not cast to ICartesianMesh*");
  return adud->data();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
