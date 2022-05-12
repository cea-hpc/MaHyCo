#ifndef CARTESIAN_CART_NODE_DIRECTION_MNG_H
#define CARTESIAN_CART_NODE_DIRECTION_MNG_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemGroup.h"
#include "cartesian/CartesianGridT.h"
#include "cartesian/CartStencilDirItemT.h"
#include "cartesian/CartLocalIdNumberingT.h"

namespace Cartesian {
  
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet les accès aux noeuds adjacents au noeud dans une direction
 */
/*---------------------------------------------------------------------------*/
class CartNode2NodeIdStencil : public CartLocalIdNumberingT<NodeLocalId> {
 public:
  //! Type d'une numérotation cartésienne sur les identifiants locaux
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;

  CartNode2NodeIdStencil(Integer dir, const CartesianNumbering& cart_numb)
  : CartLocalIdNumberingT<NodeLocalId>(cart_numb),
  m_dir (dir),
  m_nnodesm1_dir (cart_numb.nbItemDir(dir)-1),
  m_delta_dir (cart_numb.deltaDir(dir))
  {
  }

  //! Constructeur de recopie, potentiellement sur accélérateur
  ARCCORE_HOST_DEVICE CartNode2NodeIdStencil(const CartNode2NodeIdStencil& rhs)
  : CartLocalIdNumberingT<NodeLocalId>(rhs),
  m_dir (rhs.m_dir),
  m_nnodesm1_dir (rhs.m_nnodesm1_dir),
  m_delta_dir (rhs.m_delta_dir)
  {
  }

  //! Encapsulation d'un noeud central avec NLayer noeuds autour dans la direction
  /*
   *              + ----- + ----- + ----- + ----- + ----- + ----- +     ---> dir
   * NodeLocalId -1      nm2     nm1    [nid]    np1     np2     np3
   * ilayer      -3      -2      -1       0      +1      +2      +3
   */
  template<Integer NLayer>
  ARCCORE_HOST_DEVICE auto stencilNode(NodeLocalId nid, IdxType idx) const {
    return CartStencilDirItemT<NodeLocalId,NLayer>(nid, idx[m_dir], m_nnodesm1_dir, m_delta_dir);
  }

 private:
  Integer m_dir;  //! Direction privilegiee
  LocalIdType m_nnodesm1_dir;  //! Nb de noeuds-1 dans la direction m_dir
  LocalIdType m_delta_dir;  //! -+delta pour passer d'un noeud à son voisin précédent/suivant dans la direction m_dir
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Meme interface que NodeDirectionMng
 */
/*---------------------------------------------------------------------------*/
class CartNodeDirectionMng {
 public:

  //! Type pour grille cartésienne sur les identifiants locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = CartesianGrid::CartesianNumbering;

  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 public:
  CartNodeDirectionMng(const ItemInternalPtr* internals, 
    Integer dir, const CartesianGrid &cart_grid) 
  : m_internals (internals), 
  m_dir (dir),
  m_cart_grid (cart_grid),
  m_cart_numb_node (m_cart_grid.cartNumNode()),
  m_nnodes_dir (m_cart_numb_node.nbItem3()) {

    // Pour les groupes
    for(Integer d(0) ; d < 3 ; ++d) {
      if (d == m_dir) {
        // on retire la première et la dernière couche de nodes selon m_dir
        m_inner_nodes_beg[m_dir] = 1; 
        m_inner_nodes_end[m_dir] = m_nnodes_dir[m_dir]-1;
        // on ne retient que la première couche de nodes selon m_dir
        m_prev_outer_nodes_end[m_dir] = 1; 
        // on ne retient que la dernière couche de nodes selon m_dir
        m_next_outer_nodes_beg[m_dir] = m_nnodes_dir[m_dir]-1;
      } else {
        // On prend l'intégralité du domaine selon m_dir
        m_inner_nodes_beg[d] = 0;
        m_inner_nodes_end[d] = m_nnodes_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_prev_outer_nodes_end[d] = m_nnodes_dir[d];
        // On prend l'intégralité du domaine selon m_dir
        m_next_outer_nodes_beg[d] = 0;
      }
    }
  }

  CartNodeDirectionMng(const CartNodeDirectionMng& rhs)
  : m_internals (rhs.m_internals),
  m_dir (rhs.m_dir),
  m_cart_grid (rhs.m_cart_grid),
  m_cart_numb_node (rhs.m_cart_numb_node),
  m_nnodes_dir (rhs.m_nnodes_dir)
  {
    for(Integer d(0) ; d < 3 ; ++d) {
      m_inner_nodes_beg[d] = rhs.m_inner_nodes_beg[d];
      m_inner_nodes_end[d] = rhs.m_inner_nodes_end[d];
      m_prev_outer_nodes_end[d] = rhs.m_prev_outer_nodes_end[d];
      m_next_outer_nodes_beg[d] = rhs.m_next_outer_nodes_beg[d];
    }
  }

  //! Création d'une instance de Node à partir de son local_id
  Node toNode(LocalIdType node_id) const {
    if (ItemId::null(node_id)) {
      return Node();
    } else {
      return Node(m_internals, node_id);
    }
  }

  //! Pour passer d'un noeuds à ces noeuds voisins dans la direction
  auto node2NodeIdStencil() const {
    return CartNode2NodeIdStencil(m_dir, m_cart_numb_node);
  }

  //! Retourne le groupe de toutes les noeuds cartésiens
  CartNodeGroup allNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, {0, 0, 0}, m_nnodes_dir);
  }

  //! Groupe de tous les noeuds cartesiens internes à la direction
  CartNodeGroup innerNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, m_inner_nodes_beg, m_inner_nodes_end);
  }

  //! Groupe de tous les noeuds cartesiens externes à la direction à gauche (le noeud avant est nul)
  CartNodeGroup previousOuterNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, {0, 0, 0}, m_prev_outer_nodes_end);
  }

  //! Groupe de tous les noeuds cartesiens externes à la direction à droite (le noeud après est nul)
  CartNodeGroup nextOuterNodes() const {
    return CartNodeGroup(m_internals, m_dir, m_cart_grid, m_cart_numb_node, m_next_outer_nodes_beg, m_nnodes_dir);
  }

  eMeshDirection direction() const {
    return eMeshDirection(m_dir);
  }

 private:
  const ItemInternalPtr* m_internals;  //! Tableau dimensionne au nb total de Node, chaque case pointe vers un ItemInternal
  Integer m_dir;  //! Direction privilegiee

  const CartesianGrid &m_cart_grid;  //! Grille cartésienne contenant toutes les numérotations des différents items
  const CartesianNumbering &m_cart_numb_node; //! Numérotation cartésienne des noeuds (ids locaux)

  const LocalIdType3 &m_nnodes_dir; //! Nb de noeuds par direction

  LocalIdType3 m_inner_nodes_beg;  //! Triplet inclu "en bas à gauche" délimitant les noeuds intérieurs à la direction m_dir
  LocalIdType3 m_inner_nodes_end;  //! Triplet exclu "en haut à droite" délimitant les noeuds intérieurs à la direction m_dir
  LocalIdType3 m_prev_outer_nodes_end;  //! Triplet exclu "en haut à droite" délimitant les noeuds de bord au début de la direction m_dir
  LocalIdType3 m_next_outer_nodes_beg;  //! Triplet inclu "en bas à gauche" délimitant les noeuds de bord à la fin de la direction m_dir
};

}

#endif

