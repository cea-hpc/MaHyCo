#ifndef CARTESIAN_CART_CELL_DIRECTION_MNG_H
#define CARTESIAN_CART_CELL_DIRECTION_MNG_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemGroup.h"
#include "cartesian/CartesianGridT.h"
#include "cartesian/CartStencilDirItemT.h"
#include "cartesian/CartLocalIdNumberingT.h"

namespace Cartesian {
  
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Meme interface que DirCell mais en cartésien
 */
/*---------------------------------------------------------------------------*/
class CartDirCell {
 public:
  using CellType = CellLocalId;

  ARCCORE_HOST_DEVICE CartDirCell(LocalIdType cell_id, LocalIdType idx_dir, LocalIdType ncellsm1_dir, LocalIdType delta_dir)
  : m_cell_id (cell_id), 
  m_idx_dir (idx_dir),
  m_ncellsm1_dir (ncellsm1_dir),
  m_delta_dir (delta_dir) {
  }

  ARCCORE_HOST_DEVICE CellType previous() const {
    return CellType(m_idx_dir == 0 ? -1 : m_cell_id - m_delta_dir);
  }

  ARCCORE_HOST_DEVICE CellType next() const {
    return CellType(m_idx_dir == m_ncellsm1_dir ? -1 : m_cell_id + m_delta_dir);
  }

 private:

  LocalIdType m_cell_id; // le local id de la maille
  LocalIdType m_idx_dir; // indice cartesien de la maille dans la direction dir
  LocalIdType m_ncellsm1_dir; // Nb de mailles -1 dans la direction dir
  LocalIdType m_delta_dir; // +-delta a appliquer sur m_cell_id pour passer a la maille suivante/precedente
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Meme interface que CartDirCell mais avec hypothese que la maille cartésienne est interne
 * (les mailles previous()/next()) existent forcement
 */
/*---------------------------------------------------------------------------*/
class InnerCartDirCell {
 public:
  using CellType = CellLocalId;

  InnerCartDirCell(LocalIdType cell_id, LocalIdType delta_dir)
  : m_cell_id (cell_id), 
  m_delta_dir (delta_dir) {
  }

  CellType previous() const {
    return CellType(m_cell_id - m_delta_dir);
  }

  CellType next() const {
    return CellType(m_cell_id + m_delta_dir);
  }

 private:

  LocalIdType m_cell_id; // le local id de la maille
  LocalIdType m_delta_dir; // +-delta a appliquer sur m_cell_id pour passer a la maille suivante/precedente
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet les accès aux mailles adjacentes à la maille dans une direction
 */
/*---------------------------------------------------------------------------*/
class CartCell2CellIdStencil : public CartLocalIdNumberingT<CellLocalId> {
 public:
  //! Type d'une numérotation cartésienne sur les identifiants locaux
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;

  CartCell2CellIdStencil(Integer dir, const CartesianNumbering& cart_numb)
  : CartLocalIdNumberingT<CellLocalId>(cart_numb),
  m_dir (dir),
  m_ncellsm1_dir (cart_numb.nbItemDir(dir)-1),
  m_delta_dir (cart_numb.deltaDir(dir))
  {
  }

  //! Constructeur de recopie, potentiellement sur accélérateur
  ARCCORE_HOST_DEVICE CartCell2CellIdStencil(const CartCell2CellIdStencil& rhs)
  : CartLocalIdNumberingT<CellLocalId>(rhs),
  m_dir (rhs.m_dir),
  m_ncellsm1_dir (rhs.m_ncellsm1_dir),
  m_delta_dir (rhs.m_delta_dir)
  {
  }

  //! Pour accéder aux mailles voisines previous/next
  /*
   *               ------- ------- ------- 
   *              |       |       |       | 
   * CellLocalId  | prev  | [cid] |  next |   ---> dir
   *              |       |       |       |
   *               ------- ------- ------- 
   */
  ARCCORE_HOST_DEVICE auto cell(CellLocalId cid, IdxType idx) const {
    return CartDirCell(cid, idx[m_dir], m_ncellsm1_dir, m_delta_dir);
  }

  //! Encapsulation d'une maille centrale avec NLayer mailles autour dans la direction
  /*
   *               ------- ------- ------- ------- ------- ------- ------- 
   *              |       |       |       |       |       |       |       | 
   * CellLocalId  |  -1   |  cm2  |  cm1  | [cid] |  cp1  |  cp2  |  cp3  |   ---> dir
   *              |       |       |       |       |       |       |       |
   *               ------- ------- ------- ------- ------- ------- ------- 
   * ilayer          -3      -2      -1       0      +1      +2      +3
   */
  template<Integer NLayer>
  ARCCORE_HOST_DEVICE auto stencilCell(CellLocalId cid, IdxType idx) const {
    return CartStencilDirItemT<CellLocalId,NLayer>(cid, idx[m_dir], m_ncellsm1_dir, m_delta_dir);
  }

 private:
  Integer m_dir;  //! Direction privilegiee
  LocalIdType m_ncellsm1_dir;  //! Nb de noeuds-1 dans la direction m_dir
  LocalIdType m_delta_dir;  //! -+delta pour passer d'un noeud à son voisin précédent/suivant dans la direction m_dir
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet les accès aux faces adjacentes à la maille dans une direction 
 */
/*---------------------------------------------------------------------------*/
class CartDirCellFace {
 public:
  using CellType = CellLocalId;
  using FaceType = FaceLocalId;

  ARCCORE_HOST_DEVICE CartDirCellFace(LocalIdType cell_id, 
      LocalIdType prev_face_id, LocalIdType delta_next_face) 
  : m_cell_id (cell_id), 
  m_prev_face_id(prev_face_id), 
  m_delta_next_face(delta_next_face) {
  }

  ARCCORE_HOST_DEVICE CellType cellId() const {
    return CellType{m_cell_id};
  }

  ARCCORE_HOST_DEVICE FaceType previousId() const {
    return FaceType{m_prev_face_id};
  }

  ARCCORE_HOST_DEVICE FaceType nextId() const {
    return FaceType{m_prev_face_id+m_delta_next_face};
  }

 private:

  LocalIdType m_cell_id;  //! La maille pour laquelle on cherche les faces adjacentes
  LocalIdType m_prev_face_id;  //! L'identifiant de la face précédente dans la direction
  LocalIdType m_delta_next_face;  //! Delta à ajouter pour passer de la face précédente à suivante
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet d'accéder aux faces dans un stencil autour d'une maille
 * Pour passer des mailles aux faces dans un voisinage directionnel
 */
/*---------------------------------------------------------------------------*/
class CartCell2FaceIdStencil : public CartLocalIdNumberingT<CellLocalId> {
 public:
  //! Type d'une numérotation cartésienne sur les identifiants locaux
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;

  CartCell2FaceIdStencil(Integer dir, const CartesianNumbering& cart_numb_cell,
      const CartesianNumbering& cart_numb_face)
  : CartLocalIdNumberingT<CellLocalId>(cart_numb_cell),
  m_cart_numb_face_id (cart_numb_face),
  m_dir (dir),
  m_nfacesm1_dir (cart_numb_face.nbItemDir(dir)-1),
  m_delta_dir (cart_numb_face.deltaDir(dir))
  {
  }

  //! Constructeur de recopie, potentiellement sur accélérateur
  ARCCORE_HOST_DEVICE CartCell2FaceIdStencil(const CartCell2FaceIdStencil& rhs)
  : CartLocalIdNumberingT<CellLocalId>(rhs),
  m_cart_numb_face_id (rhs.m_cart_numb_face_id),
  m_dir (rhs.m_dir),
  m_nfacesm1_dir (rhs.m_nfacesm1_dir),
  m_delta_dir (rhs.m_delta_dir)
  {
  }

  //! Items adjacents à la maille de local id cid et d'indices cartésiens cidx
  ARCCORE_HOST_DEVICE CartDirCellFace cellFace(CellLocalId cid, IdxType cidx) const {
    LocalIdType prev_fid = m_cart_numb_face_id.id(cidx[0], cidx[1], cidx[2]);
    return CartDirCellFace{cid, prev_fid, m_delta_dir};
  }

  //! Encapsulation d'une maille centrale avec NLayer faces autour dans la direction
  /*
   *               ------- ------- ------- ------- ------- 
   *              |       |       |       |       |       | 
   * FaceLocalId -1      fm2     fm1[cid]fp1     fp2     fp3   ---> dir
   *              |       |       |       |       |       |
   *               ------- ------- ------- ------- ------- 
   * ilayer      -3      -2      -1      +1      +2      +3
   *
   */
  template<Integer NLayer>
  ARCCORE_HOST_DEVICE auto stencilCell2Face(CellLocalId cid, IdxType cidx) const {
    // A partir d'une maille (i,j,k) on détermine l'id de la face avec le même (i,j,k)
    // Cette face se trouve juste avant la maille d'où PREV_fid
    LocalIdType prev_fid = m_cart_numb_face_id.id(cidx[0], cidx[1], cidx[2]); 
    return NegAsymStencilDirItemT<FaceLocalId,NLayer>(prev_fid, cidx[m_dir], m_nfacesm1_dir, m_delta_dir);
  }

 private:
  CartLocalIdNumberingT<FaceLocalId> m_cart_numb_face_id;  //! Numérotation allégée aux faces
  Integer m_dir;  //! Direction privilegiee
  LocalIdType m_nfacesm1_dir;  //! Nb de faces-1 dans la direction m_dir
  LocalIdType m_delta_dir;  //! -+delta pour passer d'une face à sa voisine précédente/suivante dans la direction m_dir
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Quasi meme interface que DirCellNode
 * TODO : a supprimer et a remplacer par CartConnectivityCellFaceNode
 */
/*---------------------------------------------------------------------------*/
class CartDirCellNode {
 public:
  using NodeType = NodeLocalId;

  CartDirCellNode(LocalIdType base_node_id, const LocalIdType4 &nodef_stride)
  : m_base_node_id(base_node_id),
  m_nodef_stride(nodef_stride) {
  }

  // Nb de noeuds sur la face orthogonale a dir
  static Integer nbNode(Integer dimension) {
    return 1 << (dimension-1);
  }

  NodeType node(Integer inode) const {
    return NodeType(m_base_node_id + m_nodef_stride[inode]);
  }

 private:

  LocalIdType m_base_node_id; // le local id du noeud le plus en bas a gauche de la maille (noeud (i,j,k) de la maille (i,j,k))
  const LocalIdType4 &m_nodef_stride; // Pour une face orthogonale a dir et pour un noeud de "base", les sauts pour trouver les autres noeuds
};


/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Meme interface que CellDirectionMng
 * TODO : factoriser le calcul des noeuds voisins avec CartConnectivityCellFaceNode
 */
/*---------------------------------------------------------------------------*/
class CartCellDirectionMng {
 public:
  using CellType = CellLocalId;
  using NodeType = NodeLocalId;
  using CellEnumeratorType = CartCellEnumerator;

  //! Type pour grille cartésienne sur les identifiants locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = CartesianGrid::CartesianNumbering;

  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 public:
  CartCellDirectionMng(const ItemInternalPtr* internals, 
    Integer dir, const CartesianGrid &cart_grid) 
  : m_internals (internals), 
  m_dir (dir),
  m_cart_grid (cart_grid),
  m_cart_numb_cell (m_cart_grid.cartNumCell()),
  m_cart_numb_node (m_cart_grid.cartNumNode()),
  m_cart_numb_face_dir (m_cart_grid.cartNumFace(m_dir)),
  m_ncells_dir (m_cart_numb_cell.nbItem3()) {

    m_ncellsm1_dir = m_ncells_dir[m_dir] - 1;
    m_delta_dir = m_cart_numb_cell.deltaDir(m_dir);
    m_delta_next_face_dir = m_cart_grid.cartNumFace(m_dir).deltaDir(m_dir);

    // Pour les groupes
    for(Integer d(0) ; d < 3 ; ++d) {
      if (d == m_dir) {
        // on retire la première et la dernière couche de mailles selon m_dir
        m_inner_cells_beg[m_dir] = 1; 
        m_inner_cells_end[m_dir] = m_ncellsm1_dir;
        // on ne retient que la première couche de mailles selon m_dir
        m_prev_outer_cells_end[m_dir] = 1; 
        // on ne retient que la dernière couche de cells selon m_dir
        m_next_outer_cells_beg[m_dir] = m_ncellsm1_dir;
      } else {
        // On prend l'intégralité du domaine selon m_dir
        m_inner_cells_beg[d] = 0;
        m_inner_cells_end[d] = m_ncells_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_prev_outer_cells_end[d] = m_ncells_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_next_outer_cells_beg[d] = 0;
      }
    }

    // Pour les noeuds
    // Pour une face orthogonale a dir et pour un noeud de "base", on determine les sauts pour trouver les autres noeuds
    Integer dim = m_cart_grid.dimension();
    Integer nb_nodes_face_dir = 1 << (dim-1); // Nb de noeuds sur une face
    LocalIdType3 prev_deca[4] = {
      {0, 0, 0}, // Le premier noeud est celui de "base"
      {0, 1, 0},
      {0, 0, 1},
      {0, 1, 1},
    };
    _stride_node(nb_nodes_face_dir, prev_deca, m_nodef_stride[MS_previous]);

    LocalIdType3 next_deca[4] = {
      {1, 0, 0}, 
      {1, 1, 0},
      {1, 0, 1},
      {1, 1, 1},
    };
    _stride_node(nb_nodes_face_dir, next_deca, m_nodef_stride[MS_next]);
  }

  CartCellDirectionMng(const CartCellDirectionMng& rhs)
  : m_internals (rhs.m_internals),
  m_dir (rhs.m_dir),
  m_cart_grid (rhs.m_cart_grid),
  m_cart_numb_cell (rhs.m_cart_numb_cell),
  m_cart_numb_node (rhs.m_cart_numb_node),
  m_cart_numb_face_dir (rhs.m_cart_numb_face_dir),
  m_ncells_dir (rhs.m_ncells_dir),
  m_ncellsm1_dir (rhs.m_ncellsm1_dir),
  m_delta_dir (rhs.m_delta_dir),
  m_delta_next_face_dir (rhs.m_delta_next_face_dir)
  {
    for(Integer d(0) ; d < 3 ; ++d) {
      m_inner_cells_beg[d] = rhs.m_inner_cells_beg[d];
      m_inner_cells_end[d] = rhs.m_inner_cells_end[d];
      m_prev_outer_cells_end[d] = rhs.m_prev_outer_cells_end[d];
      m_next_outer_cells_beg[d] = rhs.m_next_outer_cells_beg[d];
    }
    for(Integer inode = 0 ; inode < 4 ; inode++) {
      m_nodef_stride[MS_previous][inode] = rhs.m_nodef_stride[MS_previous][inode];
      m_nodef_stride[MS_next][inode] = rhs.m_nodef_stride[MS_next][inode];
    }
  }

  CartDirCell cell(const CellEnumeratorType &c) const {
    return CartDirCell(c.localId(), c.itemIdx()[m_dir], m_ncellsm1_dir, m_delta_dir);
  }

  CartDirCell operator[](const CellEnumeratorType &c) const {
    return CartDirCell(c.localId(), c.itemIdx()[m_dir], m_ncellsm1_dir, m_delta_dir);
  }

  // TODO : test
  InnerCartDirCell innerCell(const CellEnumeratorType &c) const {
    return InnerCartDirCell(c.localId(), m_delta_dir);
  }

  //! Pour passer d'une maille à ces mailles voisines dans la direction 
  auto cell2CellIdStencil() const {
    return CartCell2CellIdStencil(m_dir, m_cart_numb_cell);
  }

  //! Pour passer des mailles aux faces dans un voisinage directionnel
  auto cell2FaceIdStencil() const {
    return CartCell2FaceIdStencil(m_dir, m_cart_numb_cell, m_cart_numb_face_dir);
  }

  CartDirCellNode cellNode(const CellEnumeratorType &c, eMeshSide side) const {
    const auto &cell_ijk = c.itemIdx(); // Indices du premier noeud de la maille
    LocalIdType base_node_id(m_cart_numb_node.id(cell_ijk));
    return CartDirCellNode(base_node_id, m_nodef_stride[side]);
  }

  // Pour accéder aux faces adjacentes dans la direction privilégiée
  CartDirCellFace cellFace(const CellEnumeratorType &c) const {
    return CartDirCellFace{c.localId(), c.localIdConv(m_type_face), m_delta_next_face_dir};
  }

  //! Création d'une instance de Cell à partir de son local_id
  Cell toCell(LocalIdType cell_id) const {
    if (ItemId::null(cell_id)) {
      return Cell();
    } else {
      return Cell(m_internals, cell_id);
    }
  }

  //! Retourne le groupe de toutes les mailles cartesiennes
  CartCellGroup allCells() const {
    return CartCellGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_cell, {0, 0, 0}, m_ncells_dir);
  }

  //! Groupe de toutes les mailles cartesiennes internes à la direction
  CartCellGroup innerCells() const {
    return CartCellGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_cell, m_inner_cells_beg, m_inner_cells_end);
  }

  //! Groupe de toutes les mailles cartesiennes externes à la direction à gauche (la maille avant est nulle)
  CartCellGroup previousOuterCells() const {
    return CartCellGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_cell, {0, 0, 0}, m_prev_outer_cells_end);
  }

  //! Groupe de toutes les mailles cartesiennes externes à la direction à droite (la maille après est nulle)
  CartCellGroup nextOuterCells() const {
    return CartCellGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_cell, m_next_outer_cells_beg, m_ncells_dir);
  }

  eMeshDirection direction() const {
    return eMeshDirection(m_dir);
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
  const ItemInternalPtr* m_internals;  //! Tableau dimensionne au nb total de Cell, chaque case pointe vers un ItemInternal
  Integer m_dir;  //! Direction privilegiee

  const CartesianGrid &m_cart_grid;  //! Grille cartésienne contenant toutes les numérotations des différents items
  const CartesianNumbering &m_cart_numb_cell;  //! Permet de numeroter a partir de (i,j,k) et reciproquement
  const CartesianNumbering &m_cart_numb_node; //! Numérotation cartésienne des noeuds (ids locaux)
  const CartesianNumbering &m_cart_numb_face_dir; //! Numérotation cartésienne des faces (ids locaux) de direction m_dir

  const LocalIdType3 &m_ncells_dir; //! Nb de mailles par direction

  LocalIdType m_ncellsm1_dir; //! Nb de mailles -1 dans la direction dir
  LocalIdType m_delta_dir; //! +-delta a appliquer sur m_cell_id pour passer a la maille suivante/precedente
  LocalIdType m_delta_next_face_dir;  //! +delta à appliquer sur face_id pour passer à la face suivante dans la direction

  LocalIdType3 m_inner_cells_beg;  //! Triplet inclu "en bas à gauche" délimitant les mailles intérieures à la direction m_dir
  LocalIdType3 m_inner_cells_end;  //! Triplet exclu "en haut à droite" délimitant les mailles intérieures à la direction m_dir
  LocalIdType3 m_prev_outer_cells_end;  //! Triplet exclu "en haut à droite" délimitant les mailles de bord au début de la direction m_dir
  LocalIdType3 m_next_outer_cells_beg;  //! Triplet inclu "en bas à gauche" délimitant les mailles de bord à la fin de la direction m_dir

  // Pour une face orthogonale a dir et pour un noeud de "base", on determine les sauts pour trouver les autres noeuds
  LocalIdType4 m_nodef_stride[MS_max] = {{0, 0, 0, 0}, {0, 0, 0, 0}};

  constexpr static Face* m_type_face = nullptr;  //! Permet d'utiliser la méthode localIdConv sur les Face
};

}

#endif

