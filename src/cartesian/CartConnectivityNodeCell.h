#ifndef CARTESIAN_CART_CONNECTIVITY_NODE_CELL_H
#define CARTESIAN_CART_CONNECTIVITY_NODE_CELL_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemEnumeratorT.h"
#include "cartesian/NumberingConverterT.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Connectivité noeud => mailles
 */
/*---------------------------------------------------------------------------*/
class NodeCellConnectivity {
 public:
  using CellType = CellLocalId;
 public:

  NodeCellConnectivity(LocalIdType cell_node_id, const LocalIdType3& idx3, const LocalIdType3& ncells, const LocalIdType* cell_stride, const LocalIdType3* deca) 
  : m_cell_node_id (cell_node_id), m_ncells (ncells),
  m_cell_stride (cell_stride),
  m_deca (deca)
  {
    m_idx3[0] = idx3[0];
    m_idx3[1] = idx3[1];
    m_idx3[2] = idx3[2];
  }
  //! Retourne la icell-ième maille du noeud sélectionné
  inline CellType cell(Integer icell) const {
    return CellType(_isValidStride(icell) ? m_cell_node_id + m_cell_stride[icell] : -1);
  }

 private:
  inline bool _isValidIdx(const LocalIdType idx, const LocalIdType nidx) const {
    return (idx>=0 && idx<nidx);
  }
  inline bool _isValidStride(Integer icell) const {
    // On regarde si les indices de la maille (obtenus par décallage) restent dans la grille de mailles
    // On fait le test en 3D pour éviter des conditions, mais au prix de plus de calcul
    return _isValidIdx(m_idx3[0]+m_deca[icell][0], m_ncells[0]) &&
      _isValidIdx(m_idx3[1]+m_deca[icell][1], m_ncells[1]) &&
      _isValidIdx(m_idx3[2]+m_deca[icell][2], m_ncells[2]);
  }

 private:
  LocalIdType m_cell_node_id;  //! le local id de la maille (i,j,k) pour un noeud (i,j,k)
  LocalIdType3 m_idx3;  //! (i,j,k)
  const LocalIdType3& m_ncells;  //! Nb de mailles dans chaque direction dans la grille cartésienne
  const LocalIdType* m_cell_stride;  //! A partir de m_cell_node_id, sauts pour trouver les autres noeuds de la maille
  const LocalIdType3* m_deca;  //! Décallages par rapport à (i,j,k) à partir desquels on a obtenu m_cell_stride
};  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Connectivité noeud internes => mailles (toutes les mailles existent)
 * TODO : TEST TEST TEST
 */
/*---------------------------------------------------------------------------*/
class InnerNodeCellConnectivity {
 public:
  using CellType = CellLocalId;
 public:

  InnerNodeCellConnectivity(LocalIdType cell_node_id, const LocalIdType* cell_stride) 
  : m_cell_node_id (cell_node_id),
  m_cell_stride (cell_stride)
  {
  }
  //! Retourne la icell-ième maille du noeud sélectionné
  inline CellType cell(Integer icell) const {
    return CellType(m_cell_node_id + m_cell_stride[icell]);
  }

 private:
  LocalIdType m_cell_node_id;  //! le local id de la maille (i,j,k) pour un noeud (i,j,k)
  const LocalIdType* m_cell_stride;  //! A partir de m_cell_node_id, sauts pour trouver les autres noeuds de la maille
};  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Connectivité noeud => mailles
 */
/*---------------------------------------------------------------------------*/
class CartConnectivityNodeCell {
 public:
  using NodeEnumeratorType = CartNodeEnumerator;
  using CartesianGrid = CartesianGridT<LocalIdType>;
  using CartesianNumbering = CartesianGrid::CartesianNumbering;

  //! Type pour convertir un numéro cartésien de noeud en numéro de maille
  using NumbConv = NumberingConverterT<Node,Cell>;

  //! Type de tri pour la liste des mailles
  enum eSortCell {
    SC_cart  = 0,  //! Tri selon la numérotation cartésienne croissante
    SC_trigo = 1,  //! Classement selon le sens trigonométrique dans le plan (X,Y) pour Z=Zmin puis Z=Zmax
    SC_arc = 2,  //! Classement selon le sens détermininé par Arcane
    SC_invalid = -1
  };
  
 public:
  CartConnectivityNodeCell(const CartesianGrid &cart_grid, eSortCell sort_cell = SC_cart)
  : m_sort_cell (sort_cell),
  m_cart_grid (cart_grid),
  m_cart_numb_cell (m_cart_grid.cartNumCell()),
  m_cart_numb_node (m_cart_grid.cartNumNode()),
  m_ncells (m_cart_numb_cell.nbItem3()),
  m_numb_conv (0, cart_grid) {

    Integer dim = m_cart_grid.dimension();
    m_max_nb_cell = 1 << dim; // 2 en 1D, 4 en 2D et 8 en 3D au max, sachant que des noeuds sur le bord n'ont pas ce nb max de mailles

    // Pour un noeud (i,j,k), la maille de base est la maille (i,j,k)
    /* En 2D
     * +---------+---------+
     * |         |         |
     * | i-1,j   |  i,j    |
     * |         |         |
     * +--------i,j--------+
     * |         |         |
     * | i-1,j-1 |  i,j-1  |
     * |         |         |
     * +---------+---------+
     *
     * En 3D
     * Plan Zinf (vers les Z<0)          Plan Zsup (vers les Z>0, dit topZ)
     * +-----------+-----------+         +-----------+-----------+
     * |           |           |         |           |           |
     * |i-1,j,k-1  |  i,j,k-1  |         |  i-1,j,k  |   i,j,k   |
     * |           |           |         |           |           |
     * +---------i,j,k---------+         +---------i,j,k---------+
     * |           |           |         |           |           |
     * |i-1,j-1,k-1| i,j-1,k-1 |         | i-1,j-1,k |  i,j-1,k  |
     * |           |           |         |           |           |
     * +-----------+-----------+         +-----------+-----------+
     *
     */
    // Pour une maille de "base", on determine les sauts pour trouver les autres mailles
    if (m_sort_cell == SC_cart && dim == 2) {
      LocalIdType3 deca[4] = {
        {-1, -1, 0},  
        { 0, -1, 0},  // 2  |  3 
        {-1,  0, 0},  // --i,j--
        { 0,  0, 0}   // 0  |  1
      };
      _cpyDeca(m_max_nb_cell, deca);

    } else if (m_sort_cell == SC_cart && dim == 3) {
      LocalIdType3 deca[8] = {
        {-1, -1, -1},   // Zinf        
        { 0, -1, -1},   // 2  |  3     
        {-1,  0, -1},   // --i,j--     
        { 0,  0, -1},   // 0  |  1     
        {-1, -1,  0},   // Zsup (topZ) 
        { 0, -1,  0},   // 6  |  7     
        {-1,  0,  0},   // --i,j--     
        { 0,  0,  0}    // 4  |  5     
      };
      _cpyDeca(m_max_nb_cell, deca);

    } else if (m_sort_cell == SC_trigo && dim == 2) {
      LocalIdType3 deca[4] = {
        {-1, -1, 0},  
        { 0, -1, 0},  // 3  |  2 
        { 0,  0, 0},  // --i,j--
        {-1,  0, 0}   // 0  |  1
      };
      _cpyDeca(m_max_nb_cell, deca);

    } else if (m_sort_cell == SC_trigo && dim == 3) {
      LocalIdType3 deca[8] = {
        {-1, -1, -1},   // Zinf        
        { 0, -1, -1},   // 3  |  2     
        { 0,  0, -1},   // --i,j--     
        {-1,  0, -1},   // 0  |  1     
        {-1, -1,  0},   // Zsup (topZ) 
        { 0, -1,  0},   // 7  |  6     
        { 0,  0,  0},   // --i,j--     
        {-1,  0,  0}    // 4  |  5     
      };
      _cpyDeca(m_max_nb_cell, deca);

    } else if (m_sort_cell == SC_arc && dim == 2) {
      LocalIdType3 deca[4] = {
        {-1,  0, 0},  
        { 0,  0, 0},  // 0  |  1 
        { 0, -1, 0},  // --i,j--
        {-1, -1, 0}   // 3  |  2
      };
      _cpyDeca(m_max_nb_cell, deca);

    } else if (m_sort_cell == SC_arc && dim == 3) {
      LocalIdType3 deca[8] = {
        {-1,  0, -1},  // Zinf       
        { 0,  0, -1},  // 0  |  1    
        { 0, -1, -1},  // --i,j--    
        {-1, -1, -1},  // 3  |  2    
        {-1,  0,  0},  // Zsup (topZ)
        { 0,  0,  0},  // 4  |  5    
        { 0, -1,  0},  // --i,j--    
        {-1, -1,  0}   // 7  |  6    
      };
      _cpyDeca(m_max_nb_cell, deca);

    } else {
      ARCANE_ASSERT(m_sort_cell != SC_invalid && dim!=1, ("Type de tri inconnu ou dimension 1 non supportée"));
    }
    _stride_cell(m_max_nb_cell, m_deca, m_cell_stride);
  }

  //! Nombre maximal de mailles auquel est connecté un noeud pour la dimension de la grille cartésienne
  inline Integer maxNbCell() const {
    return m_max_nb_cell;
  }

  //! Retourne l'objet permettant de récupérer les mailles connectés à un noeud interne identifié par son itérateur
  // TODO : preuve de concept, Utiliser plutôt un itérateur InnerNodeEnumerator ?
  inline InnerNodeCellConnectivity innerNodeConnectivity(const NodeEnumeratorType &n) const {
    return InnerNodeCellConnectivity(n.localIdConv(m_type_cell), m_cell_stride);
  }

  //! Retourne l'objet permettant de récupérer les mailles connectés à un noeud identifié par son itérateur
  inline NodeCellConnectivity nodeConnectivity(const NodeEnumeratorType &n) const {
    return NodeCellConnectivity(n.localIdConv(m_type_cell), n.itemIdx(), m_ncells, m_cell_stride, m_deca);
  }

  //! Retourne l'objet permettant de récupérer les mailles connectés à un noeud identifié par son local id
  inline NodeCellConnectivity nodeConnectivity(const LocalIdType nid) const {
    LocalIdType3 node_ijk;
    m_cart_numb_node.ijk(nid, node_ijk); // calcul des (i,j,k)
    LocalIdType delta=m_numb_conv.computeDelta(node_ijk[1], node_ijk[2]);
    return NodeCellConnectivity(nid + delta, node_ijk, m_ncells, m_cell_stride, m_deca);
  }

  //! Retourne l'objet permettant de récupérer les mailles connectés à un noeud identifié par son local id encapsulé
  inline NodeCellConnectivity nodeConnectivity(const NodeLocalId &nid) const {
    return nodeConnectivity(nid.localId());
  }

 private:
  using LocalIdType8 = LocalIdType[8];
  
 private:
  void _cpyDeca(Integer nb_cells, LocalIdType3 deca[]) {
    for(Integer icell = 0 ; icell < nb_cells ; icell++) {
      m_deca[icell][0] = deca[icell][0];
      m_deca[icell][1] = deca[icell][1];
      m_deca[icell][2] = deca[icell][2];
    }
    for(Integer icell = nb_cells ; icell < 8 ; icell++) {
      m_deca[icell][0] = 0;
      m_deca[icell][1] = 0;
      m_deca[icell][2] = 0;
    }
  }

  void _stride_cell(Integer nb_cells, 
      LocalIdType3 deca[], 
      LocalIdType8 cell_stride_rel) {

    for(Integer icell = 0 ; icell < nb_cells ; icell++) {
      // Equivalent à : deca[icell][0]*m_cell_delta_dir[0] + deca[icell][1]*m_cell_delta_dir[1] + deca[icell][2]*m_cell_delta_dir[2]
      cell_stride_rel[icell] = m_cart_numb_cell.id(deca[icell]) - m_cart_numb_cell.firstId();
    }
  }

 private:
  eSortCell m_sort_cell;  //! Choix du tri pour la liste des mailles connectées à un noeud
  const CartesianGrid &m_cart_grid;  //! Grille cartésienne sur les ids locaux
  const CartesianNumbering &m_cart_numb_cell;  //! Numérotation cartésienne des mailles (ids locaux)
  const CartesianNumbering &m_cart_numb_node;  //! Numérotation cartésienne des noeuds (ids locaux)
  const LocalIdType3& m_ncells;  //! Nb de mailles dans chaque direction de la grille cartésienne
  NumbConv m_numb_conv;  //! Permet de convertir un numéro cartésien de noeud (i,j,k) en numéro de maille (i,j,k)
  Integer m_max_nb_cell;  //! Nb maximal de mailles auquel est connecté un noeud pour la dimension de la grille cartésienne

  // Pour une maille de "base" (i,j,k), on determine les décallages (s0,s1,s2) pour trouver les autres mailles connectées à un noeud (i,j,k)
  LocalIdType3 m_deca[8];

  // m_cell_stride[icell] = saut correspondant au décallage m_deca[icell]
  LocalIdType8 m_cell_stride = {0, 0, 0, 0, 0, 0, 0, 0};

  constexpr static Cell* m_type_cell = nullptr;  //! Permet d'utiliser la méthode localIdConv sur les Cell
};

}

#endif

