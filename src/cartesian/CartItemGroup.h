#ifndef CARTESIAN_CART_ITEM_GROUP_H
#define CARTESIAN_CART_ITEM_GROUP_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemEnumeratorT.h"
#include "cartesian/Interval3T.h"
#include "arcane/utils/LoopRanges.h"

namespace Cartesian {
  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Encapsulation d'un ensemble d'items dans une grille cartesienne
 */
/*---------------------------------------------------------------------------*/
template<typename ITEM_TYPE>
class CartItemGoupT {
 public:
  //! Type d'itérateur correspondant aux items de type ITEM_TYPE
  using CartItemEnumerator = CartItemEnumeratorT<ITEM_TYPE>;

  //! Type de grille cartésienne cohérent avec CartItemGoupT
  using CartesianGrid = typename CartItemEnumerator::CartesianGrid;

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = typename CartItemEnumerator::CartesianNumbering;

  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = typename CartItemEnumerator::ItemInternalPtr;

  //! Type décrivant un ensemble 3D d'identifiants locaux
  using Interval3Type = Interval3T<LocalIdType>;

 public:
  CartItemGoupT(const ItemInternalPtr* internals, 
    Integer dir, const CartesianGrid &cart_grid, const CartesianNumbering &cart_item_numb, 
    const LocalIdType3 &beg, const LocalIdType3 &end) 
  : m_internals (internals),
  m_dir (dir),
  m_cart_grid (cart_grid),
  m_cart_item_numb (cart_item_numb),
  m_name ("CartItemGroup") {
    if (beg[0]<end[0] && beg[1]<end[1] && beg[2]<end[2]) {
      // L'ensemble n'est pas vide et l'on récupère les bornes
      for(Integer d(0) ; d < 3 ; ++d) {
        m_beg[d] = beg[d];
        m_end[d] = end[d];
      }
    } else {
      // L'ensemble est vide et l'on affecte tout à 0
      for(Integer d(0) ; d < 3 ; ++d) {
        m_beg[d] = 0;
        m_end[d] = 0;
      }
    }
  }

  CartItemGoupT(const CartItemGoupT& rhs)
  : m_internals (rhs.m_internals),
  m_dir (rhs.m_dir),
  m_cart_grid (rhs.m_cart_grid),
  m_name (rhs.m_name)
  {
    for(Integer d(0) ; d < 3 ; ++d) {
      m_beg[d] = rhs.m_beg[d];
      m_end[d] = rhs.m_end[d];
     }
  }

  //! Nombre d'éléments cartésiens dans le groupe
  Integer size() const {
    return (m_end[0]-m_beg[0])*(m_end[1]-m_beg[1])*(m_end[2]-m_beg[2]);
  }

  //! Nom du groupe
  const String &name() const {
    return m_name;
  }

  //! Construction d'un iterateur sur le groupe de mailles cartesiennes
  CartItemEnumerator enumerator() const {
    return CartItemEnumerator(m_internals, m_dir, m_cart_grid, m_cart_item_numb, m_beg, m_end);
  }

  //! Retourne l'intervalle 3D [m_beg[0], m_end[0][ x [m_beg[1], m_end[1][ x [m_beg[2], m_end[2][
  Interval3Type interval3() const {
    return Interval3Type(m_beg, m_end);
  }

  //! Construit {m_beg[2],size[2]} x {m_beg[1],size[1]} x {m_beg[0],size[0]} pour parcours RUNCOMMAND_LOOP
  auto loopRanges() const {
    // Attention, on donne un intervalle dans l'ordre k,j,i pour avoir un parcours cartésien
    return makeLoopRanges({m_beg[2],m_end[2]-m_beg[2]}, 
        {m_beg[1],m_end[1]-m_beg[1]}, 
        {m_beg[0],m_end[0]-m_beg[0]});
  }

 private:
  const ItemInternalPtr* m_internals;  //! Tableau dimensionne au nb total d'items ITEM_TYPE, chaque case pointe vers un ItemInternal
  Integer m_dir; //! Direction privilégiée
  const CartesianGrid &m_cart_grid; 
  const CartesianNumbering &m_cart_item_numb; //! Permet de numeroter les items

  // Ensemble des items
  // [m_beg[0], m_end[0][ x [m_beg[1], m_end[1][ x [m_beg[2], m_end[2][
  LocalIdType3 m_beg;
  LocalIdType3 m_end;

  String m_name;
};

//! Type définissant un groupe de mailles cartésiennes
using CartCellGroup = CartItemGoupT<Cell>;

// Macro generique pour parcourir avec un ItemEnumeratorT<Cell> (iterateur standard Arcane) ou un CartItemEnumeratorT<Cell>
// TODO : changer de nom ? la meme implementation peut servir pour n'importe quel item
#ifndef ENUMERATE_AUTO_CELL
#define ENUMERATE_AUTO_CELL(iname, group) for(auto iname{(group).enumerator()} ; iname.hasNext() ; ++iname)
#endif

//! Type définissant un groupe de faces cartésiennes
using CartFaceGroup = CartItemGoupT<Face>;

// Macro generique pour parcourir avec un ItemEnumeratorT<Face> (iterateur standard Arcane) ou un CartItemEnumeratorT<Face>
// TODO : changer de nom ? la meme implementation peut servir pour n'importe quel item
#ifndef ENUMERATE_AUTO_FACE
#define ENUMERATE_AUTO_FACE(iname, group) for(auto iname{(group).enumerator()} ; iname.hasNext() ; ++iname)
#endif

//! Type définissant un groupe de noeuds cartésiens
using CartNodeGroup = CartItemGoupT<Node>;

// Macro generique pour parcourir avec un ItemEnumeratorT<Node> (iterateur standard Arcane) ou un CartItemEnumeratorT<Node>
// TODO : changer de nom ? la meme implementation peut servir pour n'importe quel item
#ifndef ENUMERATE_AUTO_NODE
#define ENUMERATE_AUTO_NODE(iname, group) for(auto iname{(group).enumerator()} ; iname.hasNext() ; ++iname)
#endif

}

#endif

