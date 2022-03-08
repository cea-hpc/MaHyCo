#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/IParallelMng.h>
#include <arcane/MeshVariableScalarRef.h>
#include <arcane/MeshVariableArrayRef.h>
#include <arcane/VariableBuildInfo.h>

#include <thread>
#include <mpi.h>

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable globale         */ 
/* (i.e. non multi-mat) en utilisant des comms non-bloquantes dans des taches*/
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
void VarSyncMng::globalSynchronizeDevThr(MeshVariableRefT var)
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

#define USE_MPI_REQUEST

#ifdef USE_MPI_REQUEST
//#warning "USE_MPI_REQUEST"
  using RequestType = MPI_Request;
#else
  using RequestType = Parallel::Request;
#endif

  // L'échange proprement dit des valeurs de var
  Integer tag=1000;
  UniqueArray<RequestType> requests(2*m_nb_nei);
  IntegerUniqueArray msg_types(2*m_nb_nei); // nature des messages 

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception
    auto byte_buf_rcv_h = buf_rcv_h.byteBuf(inei); // le buffer de réception pour inei
#ifdef USE_MPI_REQUEST
    MPI_Irecv(byte_buf_rcv_h.data(), byte_buf_rcv_h.size(), MPI_BYTE, rank_nei, tag, 
        MPI_COMM_WORLD, &(requests[inei]));
#else
    requests[inei] = m_pm->recv(byte_buf_rcv_h, rank_nei, /*blocking=*/false);
#endif
    msg_types[inei] = inei+1; // >0 pour la réception
  }

  // La tâche à effectuer pour un voisin
  auto lbd_sender = [&](Integer inei, RunQueue& queue) {

    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    auto byte_buf_snd_d = buf_snd_d.byteBuf(inei); // le buffer d'envoi pour inei sur le DEVICE
    auto byte_buf_snd_h = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei sur l'HOTE

    // "byte_buf_snd_d <= var"
    async_pack_var2buf(owned_item_idx_pn[inei], var, byte_buf_snd_d, queue);

    // transfert buf_snd_d[inei] => buf_snd_h[inei]
    async_transfer(byte_buf_snd_h, byte_buf_snd_d, queue);

    // attendre que la copie sur l'hôte soit terminée pour envoyer les messages
    queue.barrier(); 

    // On amorce les envois
#ifdef USE_MPI_REQUEST
    MPI_Isend(byte_buf_snd_h.data(), byte_buf_snd_h.size(), MPI_BYTE, rank_nei, tag,
        MPI_COMM_WORLD, &(requests[m_nb_nei+inei]));
#else
    requests[m_nb_nei+inei] = m_pm->send(byte_buf_snd_h, rank_nei, /*blocking=*/false);
#endif
    msg_types[m_nb_nei+inei] = -inei-1; // <0 pour l'envoi
  };

//#define USE_THR_SENDER
#ifdef USE_THR_SENDER
#warning "USE_THR_SENDER"
  UniqueArray<std::thread*> thr_sender(m_nb_nei);

  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    thr_sender[inei] = 
      new std::thread(lbd_sender, inei, std::ref(m_neigh_queues->queue(inei)));
  }

  // On attend la fin de tous les threads
  for(auto thr : thr_sender) {
    thr->join();
    delete thr;
  }
#else
  // On lance en SEQUENTIEL
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    lbd_sender(inei, m_neigh_queues->queue(inei));
  }
#endif

  // Tache qui unpack les données du buffer reçues par un voisin inei
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst
  auto lbd_unpacker = [&](Integer inei, RunQueue& queue) {
    auto byte_buf_rcv_h = buf_rcv_h.byteBuf(inei); // buffer des données reçues sur l'HOTE
    auto byte_buf_rcv_d = buf_rcv_d.byteBuf(inei); // buffer des données reçues à transférer sur le DEVICE

    // transfert buf_rcv_h[inei] => buf_rcv_d[inei]
    async_transfer(byte_buf_rcv_d, byte_buf_rcv_h, queue);

    // "var <= buf_rcv_d[inei]"
    async_unpack_buf2var(ghost_item_idx_pn[inei], byte_buf_rcv_d, var, queue);

    // attendre que les copies soient terminées sur GPU
    queue.barrier(); 
  };
  UniqueArray<std::thread*> thr_unpacker;

  ARCANE_ASSERT(2*m_nb_nei==requests.size(), 
      ("Le nb de requetes n'est pas egal à 2 fois le nb de voisins"));

  UniqueArray<RequestType> requests2(2*m_nb_nei);
  IntegerUniqueArray msg_types2(2*m_nb_nei); 
  UniqueArray<bool> is_done_req(2*m_nb_nei);
#ifdef USE_MPI_REQUEST
  IntegerUniqueArray array_of_indices(2*m_nb_nei);
#endif

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
#ifdef USE_MPI_REQUEST
    Integer nb_req_done=0;
    MPI_Waitsome(pending_requests.size(), pending_requests.data(),
        &nb_req_done, array_of_indices.data(), MPI_STATUSES_IGNORE);
    IntegerArrayView done_indexes(array_of_indices.subView(0, nb_req_done));
#else
    IntegerUniqueArray done_indexes = m_pm->waitSomeRequests(pending_requests);
#endif

    for(Integer idone_req : done_indexes) {
      if (pending_types[idone_req] > 0) { // >0 signifie que c'est une requête de reception

        nb_pending_rcv--; // on une requete de reception en moins

        // On récupère l'indice du voisin
        Integer inei = pending_types[idone_req]-1;
        ARCANE_ASSERT(inei>=0 && inei<m_nb_nei, ("Mauvais indice de voisin"));

        // Maintenant qu'on a reçu le buffer pour le inei-ième voisin, 
        // on unpack les donnees dans un thread
//#define USE_THR_UNPACKER
#ifdef USE_THR_UNPACKER
#warning "USE_THR_UNPACKER"
        thr_unpacker.add(
            new std::thread(lbd_unpacker, inei, std::ref(m_neigh_queues->queue(inei)))
            );
#else
        lbd_unpacker(inei, m_neigh_queues->queue(inei));
#endif
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
#if 0
    std::ostringstream ostr;
    ostr << "P=" << m_pm->commRank() 
      << ", iter_wait_some=" << nb_iter_wait_some
      << ", nb_done=" << done_indexes.size();
    std::cout << ostr.str() << std::endl;
#endif
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
#if 0
    std::ostringstream ostr;
    ostr << "P=" << m_pm->commRank() 
      << ", WaitAll pending_requests.size()=" << pending_requests.size();
    std::cout << ostr.str() << std::endl;
#endif
#ifdef USE_MPI_REQUEST
    MPI_Waitall(pending_requests.size(), pending_requests.data(), MPI_STATUSES_IGNORE);
#else
    m_pm->waitAllRequests(pending_requests);
#endif
  }

#ifdef USE_THR_UNPACKER
#warning "USE_THR_UNPACKER : join"
  // On attend la fin de tous les threads unpackers
  for(auto thr : thr_unpacker) {
    thr->join();
    delete thr;
  }
#endif
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/

#define INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_THR(__MeshVariableRefT__) \
  template void VarSyncMng::globalSynchronizeDevThr(__MeshVariableRefT__ var)

INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_THR(VariableCellReal);
INST_VAR_SYNC_MNG_GLOBAL_SYNCHRONIZE_DEV_THR(VariableNodeReal3);

