#ifndef CARTESIAN_CART_TYPES_H
#define CARTESIAN_CART_TYPES_H

#include "arcane/utils/ArcaneGlobal.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Type pour les triplets cartésiens (i,j,k) et les triplets des dimensions (ni,nj,nk)
 */
/*---------------------------------------------------------------------------*/
using LocalIdType3  = Arcane::LocalIdType[3];
using UniqueIdType3 = Arcane::UniqueIdType[3];

using LocalIdType4 = Arcane::LocalIdType[4];

using IdxType = std::array<Arcane::Int64, 3>;

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Type de cote (previous ou next) pour une direction donnee
 */
/*---------------------------------------------------------------------------*/
enum eMeshSide
{
  //! Côté précédent
  MS_previous = 0,
  //! Côté suivant
  MS_next = 1,
  //! Nb maximal de côtés valides
  MS_max = 2,
  //! Côté invalide ou non initialisé
  MS_invalid = (-1)
};

}

#endif

