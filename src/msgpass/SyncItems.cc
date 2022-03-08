#include "msgpass/SyncItems.h"

#include <arcane/IVariableSynchronizer.h>
#include <arcane/IItemFamily.h>
#include <arcane/ItemGroup.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>
#include <arcane/VariableBuildInfo.h>

// Retourne le eItemKind en fonction de ItemType
template<typename ItemType>
eItemKind get_item_kind() {
  return IK_Unknown;
}

template<>
eItemKind get_item_kind<Cell>() {
  return IK_Cell;
}

template<>
eItemKind get_item_kind<Node>() {
  return IK_Node;
}

// Retourne le groupe de tous les items d'un ItemType donné
template<typename ItemType>
ItemGroupT<ItemType> get_all_items(IMesh* mesh) {
  ARCANE_ASSERT(false, ("get_all_items non implémenté pour cet ItemType"));
  return ItemGroupT<ItemType>();
}

template<>
ItemGroupT<Cell> get_all_items(IMesh* mesh) {
  return mesh->allCells();
}

template<>
ItemGroupT<Node> get_all_items(IMesh* mesh) {
  return mesh->allNodes();
}

// Retourne le groupe des items "own" d'un ItemType donné
template<typename ItemType>
ItemGroupT<ItemType> get_own_items(IMesh* mesh) {
  ARCANE_ASSERT(false, ("get_own_items non implémenté pour cet ItemType"));
  return ItemGroupT<ItemType>();
}

template<>
ItemGroupT<Cell> get_own_items(IMesh* mesh) {
  return mesh->ownCells();
}

template<>
ItemGroupT<Node> get_own_items(IMesh* mesh) {
  return mesh->ownNodes();
}

// Retourne le nom associé ItemType donné
template<typename ItemType>
const char* get_string_items() {
  ARCANE_ASSERT(false, ("get_string_items non implémenté pour cet ItemType"));
  return "UNDEF";
}

template<>
const char* get_string_items<Cell>() {
  return "Cell";
}

template<>
const char* get_string_items<Node>() {
  return "Node";
}

/*---------------------------------------------------------------------------*/
/* Encapsule la liste des items à envoyer/recevoir pour un type d'item donné */
/*---------------------------------------------------------------------------*/
template<typename ItemType>
SyncItems<ItemType>::SyncItems(IMesh* mesh, Int32ConstArrayView neigh_ranks,
    AccMemAdviser* acc_mem_adv) :
  m_buf_owned_item_idx    (platform::getAcceleratorHostMemoryAllocator()),
  m_indexes_owned_item_pn (platform::getAcceleratorHostMemoryAllocator()),
  m_nb_owned_item_pn      (platform::getAcceleratorHostMemoryAllocator()),
  m_buf_ghost_item_idx    (platform::getAcceleratorHostMemoryAllocator()),
  m_indexes_ghost_item_pn (platform::getAcceleratorHostMemoryAllocator()),
  m_nb_ghost_item_pn      (platform::getAcceleratorHostMemoryAllocator())
{
  eItemKind item_kind = get_item_kind<ItemType>();
  IItemFamily* item_family = mesh->itemFamily(item_kind);
  IVariableSynchronizer* var_sync = item_family->allItemsSynchronizer();
  
  Integer nb_nei = neigh_ranks.size();

  // "shared" ou "owned" : les items intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les items fantômes pour lesquels on va recevoir des informations
  m_indexes_owned_item_pn.resize(nb_nei);
  m_indexes_ghost_item_pn.resize(nb_nei);
  m_nb_owned_item_pn.resize(nb_nei);
  m_nb_ghost_item_pn.resize(nb_nei);

  Integer accu_nb_owned=0;
  Integer accu_nb_ghost=0;
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    m_nb_owned_item_pn[inei] = var_sync->sharedItems(inei).size();
    m_nb_ghost_item_pn[inei] = var_sync->ghostItems(inei).size();

    m_indexes_owned_item_pn[inei] = accu_nb_owned;
    m_indexes_ghost_item_pn[inei] = accu_nb_ghost;
    
    accu_nb_owned += m_nb_owned_item_pn[inei];
    accu_nb_ghost += m_nb_ghost_item_pn[inei];
  }
  m_buf_owned_item_idx.resize(accu_nb_owned);
  m_buf_ghost_item_idx.resize(accu_nb_ghost);

  // On construit des multi-vues sur des zones allouées en mémoire managées
  MultiArray2View<Integer> owned_item_idx_pn(m_buf_owned_item_idx.view(),
      m_indexes_owned_item_pn.constView(), m_nb_owned_item_pn.constView());

  MultiArray2View<Integer> ghost_item_idx_pn(m_buf_ghost_item_idx.view(),
      m_indexes_ghost_item_pn.constView(), m_nb_ghost_item_pn.constView());

  // On va récupérer les identifiants proprement dits
  ItemGroupT<ItemType> all_items = get_all_items<ItemType>(mesh);
  GroupIndexTable& lid_to_index = *all_items.localIdToIndex().get();

  auto lids2itemidx = [&](Int32ConstArrayView lids, Int32ArrayView item_idx)
  {
    for(Integer ilid=0 ; ilid<lids.size() ; ++ilid) {
      Int32 lid = lids[ilid];
      Integer idx_in_group = lid_to_index[lid];
      ARCANE_ASSERT(idx_in_group != -1, ("idx_in_group == -1"));
      item_idx[ilid] = idx_in_group;
    }
  };

  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    lids2itemidx(var_sync->sharedItems(inei), owned_item_idx_pn[inei]);
    lids2itemidx(var_sync->ghostItems(inei) , ghost_item_idx_pn[inei]);
  }

  // Les groupes d'items
  using ItemIdType = typename ItemType::LocalIdType;
  MeshVariableScalarRefT<ItemType,Integer> item_status(VariableBuildInfo(mesh, "TemporaryItemStatus"));
  auto arr_item_status = item_status.asArray();
  // = 0 : "private" : item intérieur qui ne participe à aucune comm
  // > 0 : "shared" : item "own" dont les valeurs doivent être envoyées
  // > 0 : "ghost"  : item fantôme dont on va recevoir une valeur
  arr_item_status.fill(0);

  IntegerUniqueArray all_shared_lids;
  IntegerUniqueArray all_ghost_lids;

  auto own_items = get_own_items<ItemType>(mesh);
  all_shared_lids.reserve(own_items.size()); // on surestime
  all_ghost_lids.reserve(all_items.size()-own_items.size());

  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    auto shared_item_idx = owned_item_idx_pn[inei];
    auto ghost_item_idx  = ghost_item_idx_pn[inei];

    // Les "shared" items sont des items qui appartiennent à own() 
    // (contrairement aux ghosts)
    const Integer shared_status=inei+1; // >0
    for(Integer idx : shared_item_idx) {
      ItemIdType lid(idx);
      if (arr_item_status[lid] == 0) {
        // Ici, l'item lid n'a pas encore été traité comme un "shared" item
        arr_item_status[lid] = shared_status;
        all_shared_lids.add(idx); // cet item ne sera ajouté qu'une seule fois
      }
    }

    // Les ghosts
    const Integer ghost_status=inei-1; // <0
    for(Integer idx : ghost_item_idx) {
      ItemIdType lid(idx);
      if (arr_item_status[lid] == 0) {
        // Ici, l'item lid n'a pas encore été traité comme un "ghost" item
        arr_item_status[lid] = ghost_status;
        all_ghost_lids.add(idx); // cet item ne sera ajouté qu'une seule fois
      }
    }
  }

  IntegerUniqueArray private_item_ids;
  private_item_ids.reserve(own_items.size()); // pour minimiser le nb d'allocations dynamiques

  ENUMERATE_(ItemType, iitem, own_items) {
    if (item_status[iitem] == 0) {
      private_item_ids.add(iitem->localId());
    }
  }

  // Création des groupes
  const char* str_items = get_string_items<ItemType>();
  m_private_items= item_family->createGroup( String("Private")+str_items, private_item_ids);
  m_shared_items = item_family->createGroup( String("Shared") +str_items, all_shared_lids);
  m_ghost_items  = item_family->createGroup( String("Ghost")  +str_items, all_ghost_lids);

  ARCANE_ASSERT(own_items.size()==(m_private_items.size()+m_shared_items.size()),
      ("own != private+shared"));
  ARCANE_ASSERT((all_items.size()-own_items.size())==m_ghost_items.size(),
      ("(all-own) != ghost"));

  // "Conseil" mémoire
  acc_mem_adv->setReadMostly(m_buf_owned_item_idx   .view());
  acc_mem_adv->setReadMostly(m_indexes_owned_item_pn.view());
  acc_mem_adv->setReadMostly(m_nb_owned_item_pn     .view());
  acc_mem_adv->setReadMostly(m_buf_ghost_item_idx   .view());
  acc_mem_adv->setReadMostly(m_indexes_ghost_item_pn.view());
  acc_mem_adv->setReadMostly(m_nb_ghost_item_pn     .view());

  // "Conseil" mémoire pour les groupes d'items
  acc_mem_adv->setReadMostly(m_private_items.view().localIds());
  acc_mem_adv->setReadMostly(m_shared_items .view().localIds());
  acc_mem_adv->setReadMostly(m_ghost_items  .view().localIds());
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_SYNC_ITEMS(__ItemType__) \
  template class SyncItems<__ItemType__>

INST_SYNC_ITEMS(Cell);
INST_SYNC_ITEMS(Node);

