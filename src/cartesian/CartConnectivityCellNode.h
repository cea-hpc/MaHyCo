#ifndef CARTESIAN_CART_CONNECTIVITY_CELL_NODE_H
#define CARTESIAN_CART_CONNECTIVITY_CELL_NODE_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemEnumeratorT.h"
#include "cartesian/NumberingConverterT.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Connectivité maille => noeuds
 */
/*---------------------------------------------------------------------------*/
class CellNodeConnectivity {
 public:
  using NodeType = NodeLocalId;
 public:

  CellNodeConnectivity(LocalIdType base_node_id, const LocalIdType* node_stride) 
  : m_base_node_id (base_node_id),
  m_node_stride (node_stride)
  {
  }
  //! Retourne le inode-ième noeud de la maille sélectionnée 
  inline NodeType node(Integer inode) const {
    return NodeType(m_base_node_id + m_node_stride[inode]);
  }
 private:
  LocalIdType m_base_node_id;  //! le local id du noeud le plus en bas a gauche de la maille (noeud (i,j,k) de la maille (i,j,k)) initialisée par initCartCell
  const LocalIdType* m_node_stride;  //! A partir de m_base_node_id, sauts pour trouver les autres noeuds de la maille
};  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Connectivité maille => noeuds
 */
/*---------------------------------------------------------------------------*/
class CartConnectivityCellNode {
 public:
  using CellEnumeratorType = CartCellEnumerator;
  using CartesianGrid = CartesianGridT<LocalIdType>;
  using CartesianNumbering = CartesianGrid::CartesianNumbering;

  //! Type pour convertir un numéro cartésien de maille en numéro de noeud
  using NumbConv = NumberingConverterT<Cell,Node>;

  //! Type de tri pour la liste des noeuds
  enum eSortNode {
    SN_cart  = 0,  //! Tri selon la numérotation cartésienne croissante
    SN_trigo = 1,  //! Classement selon le sens trigonométrique dans le plan (X,Y) pour Z=Zmin puis Z=Zmax
    SN_arc = 2,  //! Classement selon le sens détermininé par Arcane
    SN_invalid = -1
  };
  
 public:
  CartConnectivityCellNode(const CartesianGrid &cart_grid, eSortNode sort_node = SN_cart)
  : m_sort_node (sort_node),
  m_cart_grid (cart_grid),
  m_cart_numb_cell (m_cart_grid.cartNumCell()),
  m_cart_numb_node (m_cart_grid.cartNumNode()),
  m_numb_conv (0, cart_grid) {

    Integer dim = m_cart_grid.dimension();
    m_nb_node = 1 << dim; //! Nombre de noeuds d'une maille quelconque pour la dimension de la grille cartésienne

    // Pour les noeuds
    // Pour un noeud de "base", on determine les sauts pour trouver les autres noeuds
    if (m_sort_node == SN_cart) {
      LocalIdType3 deca[8] = {
        {0, 0, 0}, // Le premier noeud est celui de "base"
        {1, 0, 0},
        {0, 1, 0}, 
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1}, 
        {1, 1, 1}
      };
      _stride_node(m_nb_node, deca, m_node_stride);
    } else if (m_sort_node == SN_trigo) {
      LocalIdType3 deca[8] = {
        {0, 0, 0}, // Le premier noeud est celui de "base"
        {1, 0, 0},
        {1, 1, 0}, 
        {0, 1, 0},
        {0, 0, 1}, 
        {1, 0, 1},
        {1, 1, 1}, 
        {0, 1, 1},
      };
      _stride_node(m_nb_node, deca, m_node_stride);
    } else if (m_sort_node == SN_arc) {
      LocalIdType3 deca[8] = {
        {0, 1, 0},  // upperLeft
        {1, 1, 0},  // upperRight
        {1, 0, 0},  // lowerRight
        {0, 0, 0},  // lowerLeft
        {0, 1, 1},  // TopZ upperLeft
        {1, 1, 1},  // TopZ upperRight
        {1, 0, 1},  // TopZ lowerRight
        {0, 0, 1},  // TopZ lowerLeft
      };
      _stride_node(m_nb_node, deca, m_node_stride);
    } else {
      ARCANE_ASSERT(m_sort_node != SN_invalid, (""));
    }
  }

  // Nombre de noeuds d'une maille quelconque pour la dimension de la grille cartésienne
  inline Integer nbNode() const {
    return m_nb_node;
  }

  //! Retourne l'objet permettant de récupérer les noeuds de la maille identifiée par son itérateur
  inline CellNodeConnectivity cellConnectivity(const CellEnumeratorType &c) const {
    return CellNodeConnectivity(c.localIdConv(m_type_node), m_node_stride);
  }

  //! Retourne l'objet permettant de récupérer les noeuds de la maille identifiée par son local id
  inline CellNodeConnectivity cellConnectivity(const LocalIdType cid) const {
    LocalIdType3 cell_ijk;
    m_cart_numb_cell.ijk(cid, cell_ijk); // calcul des (i,j,k)
    LocalIdType delta=m_numb_conv.computeDelta(cell_ijk[1], cell_ijk[2]);
    return CellNodeConnectivity(cid + delta, m_node_stride);
  }

  //! Retourne l'objet permettant de récupérer les noeuds de la maille identifiée par son local id encapsulé
  inline CellNodeConnectivity cellConnectivity(const CellLocalId &cid) const {
    return cellConnectivity(cid.localId());
  }

 private:
  using LocalIdType8 = LocalIdType[8];
  
 private:
  void _stride_node(Integer nb_nodes, 
      LocalIdType3 deca[8], 
      LocalIdType8 node_stride_rel) {

    for(Integer inode = 0 ; inode < nb_nodes ; inode++) {
      // Equivalent à : deca[inode][0]*m_node_delta_dir[0] + deca[inode][1]*m_node_delta_dir[1] + deca[inode][2]*m_node_delta_dir[2]
      node_stride_rel[inode] = m_cart_numb_node.id(deca[inode]) - m_cart_numb_node.firstId();
    }
  }

 private:
  eSortNode m_sort_node;  //! Choix du tri pour la liste des noeuds d'une maille
  const CartesianGrid &m_cart_grid;  //! Grille cartésienne sur les ids locaux
  const CartesianNumbering &m_cart_numb_cell;  //! Numérotation cartésienne des mailles (ids locaux)
  const CartesianNumbering &m_cart_numb_node;  //! Numérotation cartésienne des noeuds (ids locaux)
  NumbConv m_numb_conv;  //! Permet de convertir un numéro cartésien de maille en son numéro de premier noeud
  Integer m_nb_node;  //! Nb de noeuds d'une maille quelconque pour la dimension de la grille cartésienne

  // Pour un noeud de "base", on determine les sauts pour trouver les autres noeuds d'une maille
  LocalIdType8 m_node_stride = {0, 0, 0, 0, 0, 0, 0, 0};

  constexpr static Node* m_type_node = nullptr;  //! Permet d'utiliser la méthode localIdConv sur les Node
};

}

#endif

