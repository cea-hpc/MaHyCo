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
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
Ref<GlobalSyncRequest<ItemType, DataType, MeshVarRefT> > VarSyncMng::iGlobalSynchronizeQueue(
    Ref<RunQueue> ref_queue, MeshVarRefT<ItemType, DataType> var)
{
  SyncItems<ItemType>* sync_items = getSyncItems<ItemType>();

  return makeRef(new GlobalSyncRequest<ItemType, DataType, MeshVarRefT>(sync_items, m_sync_buffers,
      m_pm, m_neigh_ranks,
      ref_queue, var));
}

template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
GlobalSyncRequest<ItemType, DataType, MeshVarRefT>::GlobalSyncRequest(
    SyncItems<ItemType>* sync_items,
    SyncBuffers* sync_buffers,
    IParallelMng* pm,
    Int32ConstArrayView neigh_ranks,
    Ref<RunQueue> ref_queue,
    MeshVarRefT<ItemType, DataType> var) :
  m_is_over     (false),
  m_pm          (pm),
  m_neigh_ranks (neigh_ranks),
  m_ref_queue   (ref_queue),
  m_var         (var)
{
  auto nb_owned_item_idx_pn = sync_items->nbOwnedItemIdxPn();
  auto nb_ghost_item_idx_pn = sync_items->nbGhostItemIdxPn();

  auto owned_item_idx_pn = sync_items->ownedItemIdxPn();
  m_ghost_item_idx_pn = sync_items->ghostItemIdxPn();

  // Pour un ItemType donné, combien de DataType sont utilisés ? => degree
  Integer degree = get_var_degree(var);

  sync_buffers->resetBuf();
  // On prévoit une taille max du buffer qui va contenir tous les messages
  sync_buffers->addEstimatedMaxSz<DataType>(nb_owned_item_idx_pn, degree);
  sync_buffers->addEstimatedMaxSz<DataType>(nb_ghost_item_idx_pn, degree);
  // Le buffer de tous les messages est réalloué si pas assez de place
  sync_buffers->allocIfNeeded();

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur l'HOTE (_h et LM_HostMem)
  m_buf_snd_h = sync_buffers->multiBufView<DataType>(nb_owned_item_idx_pn, degree, 0);
  m_buf_rcv_h = sync_buffers->multiBufView<DataType>(nb_ghost_item_idx_pn, degree, 0);

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur le DEVICE (_d et LM_DevMem)
  m_buf_snd_d = sync_buffers->multiBufView<DataType>(nb_owned_item_idx_pn, degree, 1);
  m_buf_rcv_d = sync_buffers->multiBufView<DataType>(nb_ghost_item_idx_pn, degree, 1);

  // L'échange proprement dit des valeurs de var
  Integer nb_nei = m_neigh_ranks.size();
  m_requests.reserve(2*nb_nei);

  // On amorce les réceptions
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception sur l'HOTE
    auto byte_buf_rcv = m_buf_rcv_h.byteBuf(inei); // le buffer de réception pour inei
    auto req_rcv = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    m_requests.add(req_rcv);
  }

  // On enchaine sur le device : 
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst

  // On remplit les buffers sur le DEVICE
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    auto buf_snd_inei = m_buf_snd_d.byteBuf(inei); // buffer dans lequel on va écrire
    // "buf_snd[inei] <= var"
    async_pack_var2buf(owned_item_idx_pn[inei], var, buf_snd_inei, *(ref_queue.get()));
  }

  // transfert buf_snd_d => buf_snd_h
  async_transfer(m_buf_snd_h, m_buf_snd_d, *(ref_queue.get()));
}

template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
GlobalSyncRequest<ItemType, DataType, MeshVarRefT>::~GlobalSyncRequest() {
  if (!m_is_over) {
    wait();
  }
}

template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void GlobalSyncRequest<ItemType, DataType, MeshVarRefT>::wait()
{
  Integer nb_nei = m_neigh_ranks.size();
  m_ref_queue->barrier(); // attendre que la copie sur l'hôte soit terminée pour envoyer les messages

  // On amorce les envois sur l'HOTE
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce l'envoi
    auto byte_buf_snd = m_buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei
    auto req_snd = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    m_requests.add(req_snd);
  }

  m_pm->waitAllRequests(m_requests);

  // transfert buf_rcv_h => buf_rcv_d
  async_transfer(m_buf_rcv_d, m_buf_rcv_h, *(m_ref_queue.get()));

  // On recopie les valeurs reçues dans les buffers dans m_var sur le DEVICE
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    auto buf_rcv_inei = m_buf_rcv_d.byteBuf(inei); // buffer duquel on va lire les données
    // "m_var <= buf_rcv_d[inei]"
    async_unpack_buf2var(m_ghost_item_idx_pn[inei], buf_rcv_inei, m_var, *(m_ref_queue.get()));
  }
  m_ref_queue->barrier(); // attendre que les copies soient terminées sur GPU

  m_is_over = true;
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(__ItemType__, __DataType__, __MeshVarRefT__) \
  template Ref<GlobalSyncRequest<__ItemType__, __DataType__, __MeshVarRefT__> > \
   VarSyncMng::iGlobalSynchronizeQueue(Ref<RunQueue> ref_queue, __MeshVarRefT__<__ItemType__, __DataType__> var)

INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(Cell, Real, MeshVariableScalarRefT);
INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(Cell, Real, MeshVariableArrayRefT);
INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(Cell, Integer, MeshVariableScalarRefT);

INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(Node, Real, MeshVariableScalarRefT);
INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(Node, Real3, MeshVariableScalarRefT);
INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(Node, Real, MeshVariableArrayRefT);

INST_VAR_SYNC_MNG_I_GLOBAL_SYNCHRONIZE_QUEUE(Face, Real, MeshVariableArrayRefT);
