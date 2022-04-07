#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>
#include <arcane/VariableBuildInfo.h>


/*---------------------------------------------------------------------------*/
/* Equivalent à un globalSynchronizeQueueEvent mais où les adresses          */
/* utilisées par les comms sont directement celles du device (GPU-aware)     */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronizeQueueEventD(Ref<RunQueue> ref_queue, MeshVariableRefT var)
{
  if (m_nb_nei==0) {
    ref_queue->barrier();
    return;
  }
  ARCANE_ASSERT(isDeviceAware(), 
      ("Impossibilité d'utiliser des adresses dans le DEVICE pour effectuer les comms"));

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
  // sur le DEVICE (_d et LM_DevMem)
  auto buf_snd_d = m_sync_buffers->multiBufView<DataType>(nb_owned_item_idx_pn, degree, 1);
  auto buf_rcv_d = m_sync_buffers->multiBufView<DataType>(nb_ghost_item_idx_pn, degree, 1);

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests(2*m_nb_nei);
  IntegerUniqueArray msg_types(2*m_nb_nei); // nature des messages 

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception avec les adresses sur le DEVICE (GPU-aware)
    auto byte_buf_rcv = buf_rcv_d.byteBuf(inei); // le buffer de réception pour inei
    requests[inei] = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    msg_types[inei] = inei+1; // >0 pour la réception
  }

  // On enchaine sur le device : 
  //    copie de var_dev dans buf_dev 

  // On remplit les buffers sur le DEVICE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    auto buf_snd_inei = buf_snd_d.byteBuf(inei); // buffer dans lequel on va écrire
    // "buf_snd[inei] <= var"
    async_pack_var2buf(owned_item_idx_pn[inei], var, buf_snd_inei, *(ref_queue.get()));

    // On enregistre un événement pour la fin de packing pour le voisin inei
    ref_queue->recordEvent(m_pack_events[inei]);
  }

  // Maintenant qu'on a lancé de façon asynchrones tous les packing pour tous les voisins
  // on va attendre voisin après vois la terminaison des packing pour amorcer les envois
  // avec les adresses des buffers sur le DEVICE (GPU-aware)
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // Attente de la fin du transfert
    m_pack_events[inei]->wait();

    // On amorce l'envoi
    auto byte_buf_snd = buf_snd_d.byteBuf(inei); // le buffer d'envoi pour inei
    requests[m_nb_nei+inei] = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    msg_types[m_nb_nei+inei] = -inei-1; // <0 pour l'envoi
  }

  // dans les faits, le stream est déjà synchronisé mais il faut terminer la RunQueue Arcane
  ref_queue->barrier();

  // Maitenant que toutes les requêtes de comms sont amorcées, il faut les terminer
  ARCANE_ASSERT(2*m_nb_nei==requests.size(), 
      ("Le nb de requetes n'est pas egal à 2 fois le nb de voisins"));

  UniqueArray<Parallel::Request> requests2(2*m_nb_nei);
  IntegerUniqueArray msg_types2(2*m_nb_nei); 
  UniqueArray<bool> is_done_req(2*m_nb_nei);

  // On utilise des vues pour éviter de réallouer en permanence des tableaux
  ArrayView<Parallel::Request> pending_requests(requests.view());
  ArrayView<Integer> pending_types(msg_types.view());

  ArrayView<Parallel::Request> upd_pending_requests(requests2.view());
  ArrayView<Integer> upd_pending_types(msg_types2.view());

  ArrayView<Parallel::Request> tmp_pending_requests;
  ArrayView<Integer> tmp_pending_types;

  Integer nb_iter_wait_some = 0;
  Integer nb_pending_rcv = m_nb_nei;

  while(nb_pending_rcv>0) {

    Integer nb_pending_req = pending_requests.size();

    // On dimenensionne is_done_requests au nb de requêtes d'avant waitSomeRequests
    // et on initialise à false
    ArrayView<bool> is_done_requests(is_done_req.subView(0, nb_pending_req));
    for(Integer ireq=0 ; ireq<nb_pending_req ; ++ireq) {
      is_done_requests[ireq]=false;
    }

    // Attente de quelques requetes
    IntegerUniqueArray done_indexes = m_pm->waitSomeRequests(pending_requests);

    for(Integer idone_req : done_indexes) {
      if (pending_types[idone_req] > 0) { // >0 signifie que c'est une requête de reception

        nb_pending_rcv--; // on une requete de reception en moins

        // On récupère l'indice du voisin
        Integer inei = pending_types[idone_req]-1;
        ARCANE_ASSERT(inei>=0 && inei<m_nb_nei, ("Mauvais indice de voisin"));

        // Maintenant qu'on a reçu le buffer pour le inei-ième voisin, 
        // on unpacke les donnees de façon asynchrone pour le inei-ième voisin
        {
          auto byte_buf_rcv_d = buf_rcv_d.byteBuf(inei); // buffer des données reçues à unpacker sur le DEVICE

          // Maintenant que byte_buf_rcv_d est sur DEVICE on peut enclencher le unpacking des données
          // "var <= buf_rcv_d[inei]"
          async_unpack_buf2var(ghost_item_idx_pn[inei], byte_buf_rcv_d, var, *(ref_queue.get()));
        }
      }
      is_done_requests[idone_req] = true;
    }
    // Il faut créer le nouveau tableau de requêtes pending dans upd_*
    Integer upd_nb_pending_req=0;
    for(Integer ireq=0 ; ireq<nb_pending_req ; ++ireq) {
      if (!is_done_requests[ireq]) {
        upd_pending_requests[upd_nb_pending_req]=pending_requests[ireq];
        upd_pending_types   [upd_nb_pending_req]=pending_types[ireq];
        upd_nb_pending_req++;
      }
    }

    // On échange les vues pour qu'à l'itération suivante 
    // pending_requests pointe vers upd_pending_types
    tmp_pending_requests = upd_pending_requests.subView(0, upd_nb_pending_req);
    upd_pending_requests = pending_requests;
    pending_requests = tmp_pending_requests;

    tmp_pending_types = upd_pending_types.subView(0, upd_nb_pending_req);
    upd_pending_types = pending_types;
    pending_types = tmp_pending_types;

    nb_iter_wait_some++;
  }

  // Ici, toutes les requetes de receptions sont forcement terminées 
  // (condition de la boucle while précédente)
  // Mais il peut rester encore des requetes d'envoi en cours
  if (pending_requests.size()) {
    // Normalement, il ne reste que des requêtes d'envois
    ARCANE_ASSERT(pending_requests.size()<=m_nb_nei, 
        ("Il ne peut pas rester un nb de requetes d'envoi supérieur au nb de voisins"));
    for(Integer msg_type : pending_types) {
      ARCANE_ASSERT(msg_type<0, 
          ("Un message d'envoi doit avoir un type négatif ce qui n'est pas le cas"));
    }
    m_pm->waitAllRequests(pending_requests);
  }

  // on attend la terminaison de tous les unpacks asynchrones
  ref_queue->barrier();
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_QUEUE_EVENT_D(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronizeQueueEventD(Ref<RunQueue> ref_queue, __MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_QUEUE_EVENT_D(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_QUEUE_EVENT_D(VariableNodeReal3);

