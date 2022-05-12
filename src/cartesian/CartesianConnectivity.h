// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* CartesianConnectivity.h                                     (C) 2000-2014 */
/*                                                                           */
/* Informations de connectivité d'un maillage cartésien.                     */
/*---------------------------------------------------------------------------*/
#ifndef CARTESIAN_CARTESIANCONNECTIVITY_H
#define CARTESIAN_CARTESIANCONNECTIVITY_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/TraceAccessor.h"
#include "arcane/Item.h"
#include "arcane/VariableTypedef.h"
#include "cartesian/CartesianGlobal.h"

#include "cartesian/CartesianGridT.h"
#include "cartesian/CartConnectivityCellNode.h"
#include "cartesian/CartConnectivityNodeCell.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Informations de connectivité d'un maillage cartésien.
 *
 * Pour l'instant, cela ne fonctionne que pour les maillages 2D.
 *
 * Comme tous les objets liés au maillage cartésien, ces instances ne
 * sont valides que tant que la topologie du maillage n'évolue pas.
 *
 * Cette classe sert à la fois pour les connectivités 2D et les connectivités
 * 3D. Les méthodes qui commencent par topZ ne sont valides que en 3D.
 *
 * Le nom des méthodes suit la nomenclature suivante:
 * - topZ/.: pour la direction Z
 * - upper/lower: pour la direction Y
 * - left/right: pour la direction X
 *
 * Pour la connectivité des noeuds autour d'une maille de coordonnées (X0,Y0,Z0),
 * le noeud de coordonnées (X,Y,Z) se récupère comme suit:
 * - En 3D, topZ si Z>Z0, sinon pas de préfixe. en 2D, jamais de préfixe.
 * - upper si Y>Y0, lower sinon,
 * - right si X>X0, left sinon,
 *
 * Donc par exemple, si Z>Z0, Y<Y0 et X>X0, le nom de la méthode est topZLowerRight().
 * Si Z<Z0, Y>Y0 et X>X0, le nom est upperRight().
 *
 * Le fonctionnement est le même pour les connectivités des mailles autour d'un noeud.
 */
class ARCANE_CEA_EXPORT CartesianConnectivity
{
  /*!
   * \brief Type énuméré indiquant la position.
   * \warning Les valeurs exactes ne doivent pas être utilisées car elles sont
   * susceptibles de changer.
   */
  enum ePosition
  {
    P_UpperLeft = 0,
    P_UpperRight = 1,
    P_LowerRight = 2,
    P_LowerLeft = 3,

    P_TopZUpperLeft = 4,
    P_TopZUpperRight = 5,
    P_TopZLowerRight = 6,
    P_TopZLowerLeft = 7
  };

 public:
  
  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

  //! Grille cartésienne en numérotation locale
  using CartesianGrid = CartesianGridT<LocalIdType>;

 public:

  CartesianConnectivity(const ItemInternalPtr* cell_internals,const ItemInternalPtr* node_internals,const CartesianGrid& cart_grid) 
  : m_cell_internals (cell_internals), m_node_internals (node_internals),
  m_conn_nc (cart_grid, CartConnectivityNodeCell::SC_arc),
  m_conn_cn (cart_grid, CartConnectivityCellNode::SN_arc) {
  }

  //! Maille en haut à gauche du noeud \a n
  Cell upperLeft(Node n) const { return _nodeToCell(n.localId(),P_UpperLeft); }
  //! Maille en haut à droite du noeud \a n
  Cell upperRight(Node n) const { return _nodeToCell(n.localId(),P_UpperRight); }
  //! Maille en bas à droite du noeud \a n
  Cell lowerRight(Node n) const { return _nodeToCell(n.localId(),P_LowerRight); }
  //! Maille en bas à gauche du noeud \a n
  Cell lowerLeft(Node n) const { return _nodeToCell(n.localId(),P_LowerLeft); }

  //! En 3D, maille en haut à gauche du noeud \a n
  Cell topZUpperLeft(Node n) const { return _nodeToCell(n.localId(),P_TopZUpperLeft); }
  //! En 3D, maille en haut à droite du noeud \a n
  Cell topZUpperRight(Node n) const { return _nodeToCell(n.localId(),P_TopZUpperRight); }
  //! En 3D, maille en bas à droite du noeud \a n
  Cell topZLowerRight(Node n) const { return _nodeToCell(n.localId(),P_TopZLowerRight); }
  //! En 3D, maille en bas à gauche du noeud \a n
  Cell topZLowerLeft(Node n) const { return _nodeToCell(n.localId(),P_TopZLowerLeft); }

  //! Noeud en haut à gauche de la maille \a c
  Node upperLeft(Cell c) const { return _cellToNode(c.localId(),P_UpperLeft); }
  //! Noeud en haut à droite de la maille \a c
  Node upperRight(Cell c) const { return _cellToNode(c.localId(),P_UpperRight); }
  //! Noeud en bas à droite de la maille \a c
  Node lowerRight(Cell c) const { return _cellToNode(c.localId(),P_LowerRight); }
  //! Noeud en bad à gauche de la maille \a c
  Node lowerLeft(Cell c) const { return _cellToNode(c.localId(),P_LowerLeft); }

  //! En 3D, noeud au dessus en haut à gauche de la maille \a c
  Node topZUpperLeft(Cell c) const { return _cellToNode(c.localId(),P_TopZUpperLeft); }
  //! En 3D, noeud au dessus en haut à droite de la maille \a c
  Node topZUpperRight(Cell c) const { return _cellToNode(c.localId(),P_TopZUpperRight); }
  //! En 3D, noeud au dessus en bas à droite de la maille \a c
  Node topZLowerRight(Cell c) const { return _cellToNode(c.localId(),P_TopZLowerRight); }
  //! En 3D, noeud au dessus en bas à gauche de la maille \a c
  Node topZLowerLeft(Cell c) const { return _cellToNode(c.localId(),P_TopZLowerLeft); }

  // Avec des énumérateurs cartésiens
  //! Local id de la maille en haut à gauche du noeud identifié par l'itérateur cartésien \a n
  CellLocalId upperLeft(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_UpperLeft); }
  //! Local id de la maille en haut à droite du noeud identifié par l'itérateur cartésien \a n
  CellLocalId upperRight(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_UpperRight); }
  //! Local id de la maille en bas à droite du noeud identifié par l'itérateur cartésien \a n
  CellLocalId lowerRight(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_LowerRight); }
  //! Local id de la maille en bas à gauche du noeud identifié par l'itérateur cartésien \a n
  CellLocalId lowerLeft(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_LowerLeft); }

  //! En 3D, local id de la maille en haut à gauche du noeud identifié par l'itérateur cartésien \a n
  CellLocalId topZUpperLeft(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_TopZUpperLeft); }
  //! En 3D, local id de la maille en haut à droite du noeud identifié par l'itérateur cartésien \a n
  CellLocalId topZUpperRight(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_TopZUpperRight); }
  //! En 3D, local id de la maille en bas à droite du noeud identifié par l'itérateur cartésien \a n
  CellLocalId topZLowerRight(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_TopZLowerRight); }
  //! En 3D, local id de la maille en bas à gauche du noeud identifié par l'itérateur cartésien \a n
  CellLocalId topZLowerLeft(const CartNodeEnumerator& n) const { return m_conn_nc.nodeConnectivity(n).cell(P_TopZLowerLeft); }

 private:

  //! Retourne la \a icell -ième maille du noeud repéré par son local id \a nid
  inline Cell _nodeToCell(LocalIdType nid, Integer icell) const {
    const NodeCellConnectivity&& node2cell = m_conn_nc.nodeConnectivity(nid);
    auto&& cid = node2cell.cell(icell);
    return (!ItemId::null(cid.localId()) ? Cell(m_cell_internals, cid) : Cell());
  }

  //! Retourne le \a inode -ième noeud de la maille repérée par son local id \a cid
  inline Node _cellToNode(LocalIdType cid, Integer inode) const {
    const CellNodeConnectivity&& cell2node = m_conn_cn.cellConnectivity(cid);
    auto&& nid = cell2node.node(inode);
    return Node(m_node_internals, nid);
  }

 private:

  const ItemInternalPtr* m_cell_internals=nullptr;  //! Tableau dimensionne au nb total de Cell, chaque case pointe vers un ItemInternal
  const ItemInternalPtr* m_node_internals=nullptr;  //! Tableau dimensionne au nb total de Node, chaque case pointe vers un ItemInternal
  CartConnectivityNodeCell m_conn_nc;  //! Passage local id noeud à maille en hypothèse cartésienne
  CartConnectivityCellNode m_conn_cn;  //! Passage local id maille à noeud en hypothèse cartésienne
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

