#ifndef CARTESIAN_CART_CONNECTIVITY_CELL_FACE_NODE_H
#define CARTESIAN_CART_CONNECTIVITY_CELL_FACE_NODE_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemEnumeratorT.h"

namespace Cartesian {
  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Connectivité maille => noeuds sur une face pour une direction donnée et un côté donné
 * TODO : factoriser le calcul des noeuds voisins avec CartCellDirectionMng
 */
/*---------------------------------------------------------------------------*/
class CartConnectivityCellFaceNode {
 public:
  using NodeType = NodeLocalId;
  using CellEnumeratorType = CartCellEnumerator;
  using CartesianGrid = CartesianGridT<LocalIdType>;
  using CartesianNumbering = CartesianGrid::CartesianNumbering;
  
  CartConnectivityCellFaceNode(Integer dir, const CartesianGrid &cart_grid)
  : m_dir (dir),
  m_cart_grid (cart_grid),
  m_cart_numb_node (m_cart_grid.cartNumNode()) {

    Integer dim = m_cart_grid.dimension();
    m_nb_nodes_face_dir = 1 << (dim-1); //! Nb de noeuds sur la face orthogonale a dir

    // Pour les noeuds
    // Pour une face orthogonale a dir et pour un noeud de "base", on determine les sauts pour trouver les autres noeuds
    LocalIdType3 prev_deca[4] = {
      {0, 0, 0}, // Le premier noeud est celui de "base"
      {0, 1, 0},
      {0, 0, 1},
      {0, 1, 1},
    };
    _stride_node(m_nb_nodes_face_dir, prev_deca, m_nodef_stride[MS_previous]);

    LocalIdType3 next_deca[4] = {
      {1, 0, 0}, 
      {1, 1, 0},
      {1, 0, 1},
      {1, 1, 1},
    };
    _stride_node(m_nb_nodes_face_dir, next_deca, m_nodef_stride[MS_next]);
  }

  CartConnectivityCellFaceNode(const CartConnectivityCellFaceNode& rhs)
  : m_dir (rhs.m_dir),
  m_cart_grid (rhs.m_cart_grid),
  m_cart_numb_node (rhs.m_cart_numb_node),
  m_nb_nodes_face_dir (rhs.m_nb_nodes_face_dir),
  m_base_node_id (rhs.m_base_node_id),
  m_ptr_nodef_stride (rhs.m_ptr_nodef_stride)
  {
    for(Integer inode = 0 ; inode < 4 ; inode++) {
      m_nodef_stride[MS_previous][inode] = rhs.m_nodef_stride[MS_previous][inode];
      m_nodef_stride[MS_next][inode] = rhs.m_nodef_stride[MS_next][inode];
    }
  }

  // Nb de noeuds sur la face orthogonale a dir
  Integer nbNode() const {
    return m_nb_nodes_face_dir;
  }

  void initCartCell(const CellEnumeratorType &c) {
    //const auto &cell_ijk = c.itemIdx(); // Indices du premier noeud de la maille
    //auto test_base_node_id = m_cart_numb_node.id(cell_ijk);
    m_base_node_id = c.localIdConv(m_type_node); // Indices du premier noeud de la maille
  }

  void initSide(eMeshSide side) {
    m_ptr_nodef_stride = m_nodef_stride[side];
  }

  NodeType node(Integer inode) const {
    return NodeType(m_base_node_id + m_ptr_nodef_stride[inode]);
  }

 private:
  
  void _stride_node(Integer nb_nodes_face_dir, 
      LocalIdType3 rel_deca[4], 
      LocalIdType4 nodef_stride_rel) {

    Integer dir_perp_0(1), dir_perp_1(2);
    switch (m_dir) {
      case 0: dir_perp_0 = 1; dir_perp_1 = 2; break;
      case 1: dir_perp_0 = 0; dir_perp_1 = 2; break;
      case 2: dir_perp_0 = 0; dir_perp_1 = 1; break;
    }
    LocalIdType3 deca = {0, 0, 0};
    for(Integer inode = 0 ; inode < nb_nodes_face_dir ; inode++) {

      deca[m_dir] = rel_deca[inode][0]; 
      deca[dir_perp_0] = rel_deca[inode][1];
      deca[dir_perp_1] = rel_deca[inode][2];
      // Equivalent à : deca[0]*m_node_delta_dir[0] + deca[1]*m_node_delta_dir[1] + deca[2]*m_node_delta_dir[2]
      nodef_stride_rel[inode] = m_cart_numb_node.id(deca) - m_cart_numb_node.firstId();
    }
  }

 private:
  Integer m_dir; //! Direction privilegiee
  const CartesianGrid &m_cart_grid; //! Grille cartésienne sur les ids locaux
  const CartesianNumbering &m_cart_numb_node; //! Numérotation cartésienne des noeuds (ids locaux)
  Integer m_nb_nodes_face_dir; //! Nb de noeuds sur la face orthogonale a dir

  LocalIdType m_base_node_id; // le local id du noeud le plus en bas a gauche de la maille (noeud (i,j,k) de la maille (i,j,k)) initialisée par initCartCell

  // Pour une face orthogonale a dir et pour un noeud de "base", on determine les sauts pour trouver les autres noeuds
  LocalIdType4 m_nodef_stride[MS_max] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
  LocalIdType *m_ptr_nodef_stride = nullptr; //! Pointe vers m_nodef_stride[MS_previous] ou m_nodef_stride[MS_next]

  constexpr static Node* m_type_node = nullptr;  //! Permet d'utiliser la méthode localIdConv sur les Node
};

}

#endif

