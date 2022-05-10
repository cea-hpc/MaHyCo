#include "msgpass/VarSyncMng.h"

#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>
#include <arcane/VariableBuildInfo.h>

/*---------------------------------------------------------------------------*/
/* pack_var2buf */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void pack_var2buf(IntegerConstArrayView item_idx, 
    const MeshVarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("pack_var2buf à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void pack_var2buf(IntegerConstArrayView item_idx, 
    const MeshVariableScalarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf) 
{
  auto var_arr = var.asArray();

  ArrayView<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    buf_vals[i] = var_arr[lid];
  }
}

// Spécialisation pour MeshVariable***Array***RefT
template<typename ItemType, typename DataType>
void pack_var2buf(IntegerConstArrayView item_idx, 
    const MeshVariableArrayRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf) 
{
  auto var_arr = var.asArray();

  Integer degree = var.arraySize();
  // Vue sur tableau 2D [nb_item][degree]
  Array2View<DataType> buf_vals(MultiBufView::valBuf2<DataType>(buf, degree));

  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    Span<const DataType> in_var_arr  (var_arr[lid]);
    Span<DataType>       out_buf_vals(buf_vals[i]);
    for(Integer j=0 ; j<degree ; ++j) {
      out_buf_vals[j] = in_var_arr[j];
    }
  }
}

/*---------------------------------------------------------------------------*/
/* unpack_buf2var */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void unpack_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVarRefT<ItemType, DataType> &var) 
{
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("unpack_buf2var à spécifialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void unpack_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableScalarRefT<ItemType, DataType> &var)
{
  auto var_arr = var.asArray();

  ConstArrayView<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    var_arr[lid] = buf_vals[i];
  }
}

// Spécialisation pour MeshVariable***Array***RefT
template<typename ItemType, typename DataType>
void unpack_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableArrayRefT<ItemType, DataType> &var)
{
  auto var_arr = var.asArray();

  Integer degree = var.arraySize();
  // Vue sur tableau 2D [nb_item][degree]
  ConstArray2View<DataType> buf_vals(MultiBufView::valBuf2<DataType>(buf, degree));


  Integer nb_item_idx = item_idx.size();

  for(Integer i=0 ; i<nb_item_idx ; ++i) {
    LocalIdType lid{item_idx[i]};
    Span<const DataType> in_buf_vals(buf_vals[i]);
    Span<DataType>       out_var_arr(var_arr[lid]);
    for(Integer j=0 ; j<degree ; ++j) {
      out_var_arr[j] = in_buf_vals[j];
    }
  }
}

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat)                                                      */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronize(MeshVariableRefT var)
{
  using ItemType = typename MeshVariableRefT::ItemType;
  using DataType = typename MeshVariableRefT::DataType;
  /*
  // Distinguer la lecture de var/construction des buffers des messages proprement dits

  // Remplir les buffers de communications à partir de var (pack ?)
  //  ==> il faut un type pour encapsuler les buffers de communications
  //  Idealement, les données pour tous les voisins devraient être dans un seul buffer contigu
  //  pour faciliter les transferts H<->D
  // Si un seul buffer, attention à avoir des données qui sont alignées (std::align ?)

  // Envoyer/recevoir les données avec MPI

  // Lire les buffers MPI reçus pour écrire les données dans var
  */

  SyncItems<ItemType>* sync_items = getSyncItems<ItemType>();

  auto nb_owned_item_idx_pn = sync_items->nbOwnedItemIdxPn();
  auto nb_ghost_item_idx_pn = sync_items->nbGhostItemIdxPn();

  auto owned_item_idx_pn = sync_items->ownedItemIdxPn();
  auto ghost_item_idx_pn = sync_items->ghostItemIdxPn();

  // Pour un ItemType donné, combien de DataType sont utilisés ? => degree
  Integer degree = get_var_degree(var);

  m_sync_buffers->resetBuf();
  // On prévoit une taille max du buffer qui va contenir tous les messages
  m_sync_buffers->addEstimatedMaxSz<DataType>(nb_owned_item_idx_pn, degree);
  m_sync_buffers->addEstimatedMaxSz<DataType>(nb_ghost_item_idx_pn, degree);
  // Le buffer de tous les messages est réalloué si pas assez de place
  m_sync_buffers->allocIfNeeded();

  // On récupère les adresses et tailles des buffers d'envoi et de réception
  auto buf_snd = m_sync_buffers->multiBufView<DataType>(nb_owned_item_idx_pn, degree, 0);
  auto buf_rcv = m_sync_buffers->multiBufView<DataType>(nb_ghost_item_idx_pn, degree, 0);

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests;

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception
    auto byte_buf_rcv = buf_rcv.byteBuf(inei); // le buffer de réception pour inei
    auto req_rcv = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    requests.add(req_rcv);
  }

  // On remplit les buffers sur CPU, TODO : sur GPU
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    auto buf_snd_inei = buf_snd.byteBuf(inei); // buffer dans lequel on va écrire
    // "buf_snd[inei] <= var"
    pack_var2buf(owned_item_idx_pn[inei], var, buf_snd_inei);
  }

  // On amorce les envois
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce l'envoi
    auto byte_buf_snd = buf_snd.byteBuf(inei); // le buffer d'envoi pour inei
    auto req_snd = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    requests.add(req_snd);
  }

  m_pm->waitAllRequests(requests);
  requests.clear();

  // On recopie les valeurs reçues dans les buffers dans var
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    auto buf_rcv_inei = buf_rcv.byteBuf(inei); // buffer duquel on va lire les données
    // "var <= buf_rcv[inei]"
    unpack_buf2var(ghost_item_idx_pn[inei], buf_rcv_inei, var);
  }
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/
#include <arcane/utils/Real3.h>

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronize(__MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(VariableNodeReal3);

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE(VariableCellArrayReal3);

