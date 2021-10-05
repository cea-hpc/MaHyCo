#ifndef CARTESIAN_CARTESIAN_ITEM_SORTER_H
#define CARTESIAN_CARTESIAN_ITEM_SORTER_H

#include "arcane/IMesh.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet de trier des items selon un ordre cartésien
 */
/*---------------------------------------------------------------------------*/
class CartesianItemSorter {
 public:
  CartesianItemSorter(Arcane::IMesh* mesh);
  virtual ~CartesianItemSorter();

  /*! \brief Trie les faces selon un ordre cartésien
   */
  void sortFaces();

 protected:
  Arcane::IMesh* m_mesh;
};

}

#endif
