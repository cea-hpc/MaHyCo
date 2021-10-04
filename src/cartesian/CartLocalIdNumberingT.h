#ifndef CARTESIAN_CART_LOCAL_ID_NUMBERING_T_H
#define CARTESIAN_CART_LOCAL_ID_NUMBERING_T_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartesianNumberingT.h"

#include "arcane/utils/ArrayExtents.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Gestion d'une numerotation cartesienne sur une grille d'items de 
 * dimension 3 
 * Permet d'associer une itération de boucle cartésienne avec un couple 
 * (local id, (i,j,k))
 */
/*---------------------------------------------------------------------------*/
template<typename ItemLocalIdType>
class CartLocalIdNumberingT {
 public:
  //! Type d'entiers pour identifiants
  using IdType = LocalIdType;

  //! Type d'une numérotation cartésienne sur les identifiants locaux
  using CartesianNumbering = CartesianNumberingT<IdType>;

  //! Type pour identifier un item (local id, (i,j,k))
  using LocalIdIdxType = std::pair<ItemLocalIdType,IdxType>;

  CartLocalIdNumberingT(const CartesianNumbering& cart_numb) 
  : m_first_item_id (cart_numb.firstId()),
  m_coef1 (cart_numb.deltaDir(1)),
  m_coef2 (cart_numb.deltaDir(2))
  {
  }

  //! Constructeur de recopie, potentiellement sur device
  ARCCORE_HOST_DEVICE CartLocalIdNumberingT(const CartLocalIdNumberingT<ItemLocalIdType>& rhs) 
  : m_first_item_id (rhs.m_first_item_id),
  m_coef1 (rhs.m_coef1),
  m_coef2 (rhs.m_coef2)
  {
  }

  //! Passage (i,j,k) => numero
  ARCCORE_HOST_DEVICE inline IdType id(IdType i, IdType j, IdType k) const {
    return m_first_item_id + i + j*m_coef1 + k*m_coef2;
  }

  //! Retourne le couple (local id, (i,j,k)) à partir d'un itéré d'une boucle directe
  ARCCORE_HOST_DEVICE LocalIdIdxType idIdx(const ArrayBoundsIndex<3>& iter) const {
    IdxType idx = {iter.id2(),iter.id1(),iter.id0()}; // on remet dans le bon ordre
    IdType cell_id = id(idx[0], idx[1], idx[2]); // on calcule l'id selon la numerotation cartesienne

    return {ItemLocalIdType(cell_id), idx};
  }

 private:
  IdType m_first_item_id;  //! item_id = m_first_item_id + numéro_cartésien(i,j,k), permet un décallage dans la numérotation

  //! Coefficients multiplicateurs pour formule cartésienne
  IdType m_coef1;
  IdType m_coef2;
};

}

#endif
