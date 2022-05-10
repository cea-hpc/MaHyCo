#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>
#include <arcane/VariableBuildInfo.h>

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat) dont les données sont présentes sur GPU              */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronizeQueue(Ref<RunQueue> ref_queue, MeshVariableRefT var)
{
  if (m_nb_nei==0) {
    ref_queue->barrier(); 
    return;
  }

  using ItemType = typename MeshVariableRefT::ItemType;
  using DataType = typename MeshVariableRefT::DataType;

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
  // sur l'HOTE (_h et LM_HostMem)
  auto buf_snd_h = m_sync_buffers->multiBufView<DataType>(nb_owned_item_idx_pn, degree, 0);
  auto buf_rcv_h = m_sync_buffers->multiBufView<DataType>(nb_ghost_item_idx_pn, degree, 0);

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur le DEVICE (_d et LM_DevMem)
  auto buf_snd_d = m_sync_buffers->multiBufView<DataType>(nb_owned_item_idx_pn, degree, 1);
  auto buf_rcv_d = m_sync_buffers->multiBufView<DataType>(nb_ghost_item_idx_pn, degree, 1);

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests;

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception sur l'HOTE
    auto byte_buf_rcv = buf_rcv_h.byteBuf(inei); // le buffer de réception pour inei
    auto req_rcv = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    requests.add(req_rcv);
  }

  // On enchaine sur le device : 
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst

  // On remplit les buffers sur le DEVICE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    auto buf_snd_inei = buf_snd_d.byteBuf(inei); // buffer dans lequel on va écrire
    // "buf_snd[inei] <= var"
    async_pack_var2buf(owned_item_idx_pn[inei], var, buf_snd_inei, *(ref_queue.get()));
  }

  // transfert buf_snd_d => buf_snd_h
  async_transfer(buf_snd_h, buf_snd_d, *(ref_queue.get()));
  ref_queue->barrier(); // attendre que la copie sur l'hôte soit terminée pour envoyer les messages

  // On amorce les envois sur l'HOTE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce l'envoi
    auto byte_buf_snd = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei
    auto req_snd = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    requests.add(req_snd);
  }

  m_pm->waitAllRequests(requests);
  requests.clear();

  // transfert buf_rcv_h => buf_rcv_d
  async_transfer(buf_rcv_d, buf_rcv_h, *(ref_queue.get()));

  // On recopie les valeurs reçues dans les buffers dans var sur le DEVICE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    auto buf_rcv_inei = buf_rcv_d.byteBuf(inei); // buffer duquel on va lire les données
    // "var <= buf_rcv_d[inei]"
    async_unpack_buf2var(ghost_item_idx_pn[inei], buf_rcv_inei, var, *(ref_queue.get()));
  }
  ref_queue->barrier(); // attendre que les copies soient terminées sur GPU
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_QUEUE(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronizeQueue(Ref<RunQueue> ref_queue, __MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_QUEUE(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_QUEUE(VariableNodeReal3);

