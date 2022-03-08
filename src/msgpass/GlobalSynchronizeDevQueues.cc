#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>
#include <arcane/VariableBuildInfo.h>
#include <arcane/accelerator/IRunQueueStream.h>

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat) en utilisant plusieurs queues                        */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronizeDevQueues(MeshVariableRefT var)
{
  if (m_nb_nei==0) {
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

  using RequestType = Parallel::Request;

  // L'échange proprement dit des valeurs de var
  UniqueArray<RequestType> requests(2*m_nb_nei);
  IntegerUniqueArray msg_types(2*m_nb_nei); // nature des messages 

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception
    auto byte_buf_rcv_h = buf_rcv_h.byteBuf(inei); // le buffer de réception pour inei
    requests[inei] = m_pm->recv(byte_buf_rcv_h, rank_nei, /*blocking=*/false);
    msg_types[inei] = inei+1; // >0 pour la réception
  }

  // De manière asynchrone sur GPU, on packe les données et on transfère sur CPU
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    RunQueue& queue = m_neigh_queues->queue(inei); // queue du voisin asynchrone

    auto byte_buf_snd_d = buf_snd_d.byteBuf(inei); // le buffer d'envoi pour inei sur le DEVICE
    auto byte_buf_snd_h = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei sur l'HOTE

    // "byte_buf_snd_d <= var"
    async_pack_var2buf(owned_item_idx_pn[inei], var, byte_buf_snd_d, queue);

    // transfert buf_snd_d[inei] => buf_snd_h[inei]
    async_transfer(byte_buf_snd_h, byte_buf_snd_d, queue);
  }

#ifdef ARCANE_COMPILING_CUDA
  // On récupère les streams CUDA
  UniqueArray<cudaStream_t*> pack_streams(m_nb_nei);
  UniqueArray<cudaStream_t*> pack_streams2(m_nb_nei);
  UniqueArray<Integer> pack_nei(m_nb_nei);
  UniqueArray<Integer> pack_nei2(m_nb_nei);
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    auto* rq = m_neigh_queues->queue(inei)._internalStream();
    pack_streams[inei] = reinterpret_cast<cudaStream_t*>(rq->_internalImpl());
    pack_nei[inei] = inei;
  }

  Integer nb_pend_stream=m_nb_nei; // Nb de streams pending, non synchronisés
  ArrayView<cudaStream_t*> pend_streams(pack_streams.view());
  ArrayView<cudaStream_t*> upd_pend_streams(pack_streams2.view());
  ArrayView<cudaStream_t*> tmp_pend_streams;

  // Pour faire la correspondance avec le numéro du voisin
  ArrayView<Integer> pend_nei(pack_nei.view());
  ArrayView<Integer> upd_pend_nei(pack_nei2.view());
  ArrayView<Integer> tmp_pend_nei;

  UniqueArray<bool> query_streams(m_nb_nei);

  Int64 nb_iter_stream=0;
  while(nb_pend_stream>0) {
    Integer done_streams=0, upd_nb_pend_streams=0;
    for(Integer ipend=0 ; ipend<nb_pend_stream ; ++ipend) {
      query_streams[ipend]=false;  // false = pending, true = synchronisé(terminé)
    }
    // Attente active => on sort forcement avec un stream synchronisé
    while(done_streams==0) {
      // On parcourt tous les streams pour savoir lesquels sont terminés
      for(Integer ipend=0 ; ipend<nb_pend_stream ; ++ipend) {
        if (cudaStreamQuery(*pend_streams[ipend])==cudaSuccess) {
          query_streams[ipend]=true;
          done_streams++;
        }
      }
    }
    for(Integer ipend=0 ; ipend<nb_pend_stream ; ++ipend) {
      if (query_streams[ipend]==true) {
        Integer inei = pend_nei[ipend];
        // Le stream est synchronisé 
        //  => le buffer d'envoi du inei-ième est transféré sur CPU
        //  => on peut amorcer l'envoi
        Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

        auto byte_buf_snd_h = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei sur l'HOTE

        requests[m_nb_nei+inei] = m_pm->send(byte_buf_snd_h, rank_nei, /*blocking=*/false);
        msg_types[m_nb_nei+inei] = -inei-1; // <0 pour l'envoi

      } else {
        // Le stream n'est pas encore synchronisé, on continue à l'attendre
        upd_pend_streams[upd_nb_pend_streams]=pend_streams[ipend];
        upd_pend_nei[upd_nb_pend_streams]=pend_nei[ipend];
        upd_nb_pend_streams++;
      }
    }
    tmp_pend_streams = upd_pend_streams.subView(0, upd_nb_pend_streams);
    upd_pend_streams = pend_streams;
    pend_streams = tmp_pend_streams;

    tmp_pend_nei = upd_pend_nei.subView(0, upd_nb_pend_streams);
    upd_pend_nei =     pend_nei;
        pend_nei = tmp_pend_nei;

    nb_pend_stream=upd_nb_pend_streams;
    nb_iter_stream++;
  }
#if 0
  std::ostringstream ostr;
  ostr << "P=" << m_pm->commRank() 
    << ", iter_wait_stream=" << nb_iter_stream;
  std::cout << ostr.str() << std::endl;
#endif
  // Les streams CUDA sont normalement terminées mais ils faut terminer les RunQueue
  m_neigh_queues->waitAllQueues(); 
#else
#warning "cudaStream désactivé"
  m_neigh_queues->waitAllQueues(); 
  // ici, tous les buffers d'envoi sont transférés sur CPU,
  // on peut amorcer les envois
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    auto byte_buf_snd_h = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei sur l'HOTE

    requests[m_nb_nei+inei] = m_pm->send(byte_buf_snd_h, rank_nei, /*blocking=*/false);
    msg_types[m_nb_nei+inei] = -inei-1; // <0 pour l'envoi
  }
#endif

  // Tache qui unpack les données du buffer reçues par un voisin inei
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst
  auto async_unpacker = [&](Integer inei, RunQueue& queue) {
    auto byte_buf_rcv_h = buf_rcv_h.byteBuf(inei); // buffer des données reçues sur l'HOTE
    auto byte_buf_rcv_d = buf_rcv_d.byteBuf(inei); // buffer des données reçues à transférer sur le DEVICE

    // transfert buf_rcv_h[inei] => buf_rcv_d[inei]
    async_transfer(byte_buf_rcv_d, byte_buf_rcv_h, queue);

    // "var <= buf_rcv_d[inei]"
    async_unpack_buf2var(ghost_item_idx_pn[inei], byte_buf_rcv_d, var, queue);
  };

  ARCANE_ASSERT(2*m_nb_nei==requests.size(), 
      ("Le nb de requetes n'est pas egal à 2 fois le nb de voisins"));

  UniqueArray<RequestType> requests2(2*m_nb_nei);
  IntegerUniqueArray msg_types2(2*m_nb_nei); 
  UniqueArray<bool> is_done_req(2*m_nb_nei);

  // On utilise des vues pour éviter de réallouer en permanence des tableaux
  ArrayView<RequestType> pending_requests(requests.view());
  ArrayView<Integer> pending_types(msg_types.view());

  ArrayView<RequestType> upd_pending_requests(requests2.view());
  ArrayView<Integer> upd_pending_types(msg_types2.view());

  ArrayView<RequestType> tmp_pending_requests;
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
        // on unpack les donnees de façon asynchrone
        async_unpacker(inei, m_neigh_queues->queue(inei));
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
  m_neigh_queues->waitAllQueues();
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_QUEUES(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronizeDevQueues(__MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_QUEUES(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_QUEUES(VariableNodeReal3);

