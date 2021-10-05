// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* NodeDirectionMng.cc                                         (C) 2000-2016 */
/*                                                                           */
/* Infos sur les noeuds d'une direction X Y ou Z d'un maillage structuré.    */
/*---------------------------------------------------------------------------*/
#ifndef CARTESIAN_NODEDIRECTIONMNG_H
#define CARTESIAN_NODEDIRECTIONMNG_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/ArcaneTypes.h"
#include "cartesian/CartesianGlobal.h"

#include "arcane/Item.h"
#include "arcane/ItemEnumerator.h"

#include "cartesian/CartTypes.h"
#include "cartesian/CartesianNumberingT.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class NodeDirectionMng;
class ICartesianMesh;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Noeud avant et après un noeud suivant une direction.
 *
 * Les instances de cette classe sont temporaires et construites via
 * NodeDirectionMng::node().
 */
class ARCANE_CEA_EXPORT DirNode
{
  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 public:
  DirNode(LocalIdType node_id, LocalIdType idx_dir, LocalIdType nnodesm1_dir, LocalIdType delta_dir, 
    const ItemInternalPtr* node_internals, ItemInternalPtr null_node_internal)
  : m_node_id (node_id), 
  m_idx_dir (idx_dir),
  m_nnodesm1_dir (nnodesm1_dir),
  m_delta_dir (delta_dir),
  m_node_internals (node_internals),
  m_null_node_internal (null_node_internal)
  {
  }
 public:
  //! Noeud avant
  Node previous() const { return Node(m_idx_dir>0              ? m_node_internals[m_node_id-m_delta_dir] : m_null_node_internal); }
  //! Noeud après
  Node next()     const { return Node(m_idx_dir<m_nnodesm1_dir ? m_node_internals[m_node_id+m_delta_dir] : m_null_node_internal); }
 private:
  LocalIdType m_node_id; // le local id du noeud
  LocalIdType m_idx_dir; // indice cartesien du cartésien dans la direction dir
  LocalIdType m_nnodesm1_dir; // Nb de noeuds -1 dans la direction dir
  LocalIdType m_delta_dir; // +-delta a appliquer sur m_node_id pour passer au noeud suivant/precedent
  const ItemInternalPtr* m_node_internals;
  ItemInternalPtr m_null_node_internal;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Infos sur les noeuds d'une direction spécifique X,Y ou Z
 * d'un maillage structuré.
 */
class ARCANE_CEA_EXPORT NodeDirectionMng
{
  friend class CartesianMesh;
  class Impl;

 private:

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;

 public:
  
  //! Créé une instance vide. L'instance n'est pas valide tant que init() n'a pas été appelé.
  NodeDirectionMng();
  NodeDirectionMng(const NodeDirectionMng& rhs);
  ~NodeDirectionMng();

  //! Noeud direction correspondant au noeud \a n
  DirNode node(Node n)
  {
    return _node(n.localId());
  }

  //! Groupe de tous les noeuds dans la direction.
  NodeGroup allNodes() const;

  /*!
   * \brief Groupe de tous les noeuds internes dans la direction.
   *
   * Un noeud est considéré comme interne si son noeud
   * avant ou après n'est pas nul.
   */
  NodeGroup innerNodes() const;

  /*!
   * \brief Groupe de tous les noeuds externes dans la direction.
   *
   * Un noeud est considéré comme externe si son noeud
   * avant ou après est nul.
   */
  NodeGroup outerNodes() const;

  //! Noeud direction correspondant au noeud \a n.
  DirNode operator[](Node n)
  {
    return _node(n.localId());
  }

  //! Noeud direction correspondant à l'itérateur du noeud \a inode.
  DirNode operator[](NodeEnumerator inode)
  {
    return _node(inode.itemLocalId());
  }

 private:
  
  //! Noeud direction correspondant au noeud de numéro local \a local_id
  DirNode _node(Int32 local_id)
  {
    LocalIdType idx_dir; // Trouver un moyen d'eviter les ifs sans utiliser des methodes virtuelles
    if (m_dir==0)
      idx_dir=m_p_cart_numb_node->idxDir0(local_id); // couteux
    else if (m_dir==1)
      idx_dir=m_p_cart_numb_node->idxDir1(local_id); // couteux
    else 
      idx_dir=m_p_cart_numb_node->idxDir2(local_id); // couteux

    return DirNode(local_id, /*idx3[m_dir]*/idx_dir, m_nnodesm1_dir, m_delta_dir, m_node_internals.data(), m_null_node_internal);
  }

 public:

  /*!
   * \internal
   * \brief Usage interne à Arcane. Calcul les entités internes et externes.
   * Suppose que init() a été appelé.
   */
  void computeInnerAndOuterItems(const ItemGroup& items);

  /*!
   * \internal
   * Initialise l'instance.
   */
  void init(ICartesianMesh* cm,eMeshDirection dir,CartesianNumbering* p_cart_numb_node);

  /*!
   * \internal
   * Détruit les ressources associées à l'instance.
   */
  void destroy();

  //! Valeur de la direction
  eMeshDirection direction() const
  {
    return m_direction;
  }

 private:

  eMeshDirection m_direction;
  Integer m_dir;  //! m_direction casté en Integer
  Impl* m_p;

  CartesianNumbering* m_p_cart_numb_node=nullptr; //! Pointeur sur la numérotation cartésienne aux noeuds
  ItemInternalList m_node_internals;  //! Tableau dimensionne au nb total de Node, chaque case pointe vers un ItemInternal
  LocalIdType m_nnodesm1_dir=0;
  LocalIdType m_delta_dir=0;

  Node m_null_node;  //! Un noeud nul (inexistant)
  ItemInternal* m_null_node_internal=nullptr;  //! m_null_node.internal()
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

