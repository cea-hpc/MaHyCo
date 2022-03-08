#ifndef MSG_PASS_SYNC_ITEMS_H
#define MSG_PASS_SYNC_ITEMS_H

#include "accenv/AcceleratorUtils.h"
#include <arcane/IMesh.h>
#include <arcane/utils/MultiArray2.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Encapsule la liste des items à envoyer/recevoir */
/*---------------------------------------------------------------------------*/
template<typename ItemType>
class SyncItems {
 public:
  // Le type de groupe d'items associé à ItemType
  using ItemGroupType = ItemGroupT<ItemType>;
 public:
  SyncItems(IMesh* mesh, Int32ConstArrayView neigh_ranks, 
      AccMemAdviser* acc_mem_adv);
  virtual ~SyncItems() {}

  auto nbOwnedItemIdxPn() const {
    return m_nb_owned_item_pn.constView();
  }

  auto nbGhostItemIdxPn() const {
    return m_nb_ghost_item_pn.constView();
  }

  auto ownedItemIdxPn() const {
    return ConstMultiArray2View<Integer>(m_buf_owned_item_idx.constView(),
        m_indexes_owned_item_pn.constView(), m_nb_owned_item_pn.constView());
  }

  auto ghostItemIdxPn() const {
    return ConstMultiArray2View<Integer>(m_buf_ghost_item_idx.constView(),
        m_indexes_ghost_item_pn.constView(), m_nb_ghost_item_pn.constView());
  }

  // Items intérieurs qui n'interviennent pas dans les comms
  auto privateItems() const {
    return m_private_items;
  }

  // Items intérieurs dont les valeurs vont être envoyées
  auto sharedItems() const {
    return m_shared_items;
  }

  // Items fantômes dont on va recevoir des valeurs
  auto ghostItems() const {
    return m_ghost_items;
  }

 protected:
  // "shared" ou "owned" : les items intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les items fantômes pour lesquels on va recevoir des informations
  // _pn : _per_neigh, info par voisin
  
  // Nb d'items par voisin et listes des indexes de ces items
  // Owned ou Shared
  IntegerUniqueArray m_buf_owned_item_idx;
  IntegerUniqueArray m_indexes_owned_item_pn;
  IntegerUniqueArray m_nb_owned_item_pn;

  // Ghost
  IntegerUniqueArray m_buf_ghost_item_idx;
  IntegerUniqueArray m_indexes_ghost_item_pn;
  IntegerUniqueArray m_nb_ghost_item_pn;

  // Les groupes d'items
  // own() = private + shared  , private \inter shared = 0
  // all() = own() + ghost
  ItemGroupType m_private_items;
  ItemGroupType m_shared_items;
  ItemGroupType m_ghost_items;
};

/* Retourne les groupes associés aux types d'items
 */
template<typename ItemType>
ItemGroupT<ItemType> get_all_items(IMesh* mesh);

template<typename ItemType>
ItemGroupT<ItemType> get_own_items(IMesh* mesh);

#endif
