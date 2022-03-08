#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/utils/TraceInfo.h>
#include <arcane/IParallelMng.h>

#include <arcane/materials/IMeshMaterialVariableSynchronizer.h>

/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable multi-env       */ 
/*---------------------------------------------------------------------------*/
template<typename DataType>
void VarSyncMng::multiMatSynchronize(Ref<RunQueue> ref_queue, 
    CellMaterialVariableScalarRef<DataType> var_menv)
{
  if (m_nb_nei==0) {
    return;
  }

  auto nb_owned_evi_pn = m_sync_evi->nbOwnedEviPn();
  auto nb_ghost_evi_pn = m_sync_evi->nbGhostEviPn();

  auto owned_evi_pn = m_sync_evi->ownedEviPn();
  auto ghost_evi_pn = m_sync_evi->ghostEviPn();

  // Pour une Cell donnée, combien de DataType sont utilisés ? => degree
  Integer degree = var_menv.globalVariable().arraySize();
  degree = (degree==0 ? 1 : degree);
#if 0
  std::cout << String::format("P{0} : var={1}, degree={2}",
      m_pm->commRank(), var_menv.name(), degree).localstr() << std::endl;
#endif

  m_sync_buffers->resetBuf();
  // On prévoit une taille max du buffer qui va contenir tous les messages
  m_sync_buffers->addEstimatedMaxSz<DataType>(nb_owned_evi_pn, degree);
  m_sync_buffers->addEstimatedMaxSz<DataType>(nb_ghost_evi_pn, degree);
  // Le buffer de tous les messages est réalloué si pas assez de place
  m_sync_buffers->allocIfNeeded();

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur l'HOTE (_h et LM_HostMem)
  auto buf_snd_h = m_sync_buffers->multiBufView<DataType>(nb_owned_evi_pn, degree, 0);
  auto buf_rcv_h = m_sync_buffers->multiBufView<DataType>(nb_ghost_evi_pn, degree, 0);

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur le DEVICE (_d et LM_DevMem)
  auto buf_snd_d = m_sync_buffers->multiBufView<DataType>(nb_owned_evi_pn, degree, 1);
  auto buf_rcv_d = m_sync_buffers->multiBufView<DataType>(nb_ghost_evi_pn, degree, 1);

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests(2*m_nb_nei);
  IntegerUniqueArray msg_types(2*m_nb_nei); // nature des messages 

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception sur l'HOTE
    auto byte_buf_rcv = buf_rcv_h.byteBuf(inei); // le buffer de réception pour inei
    requests[inei] = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    msg_types[inei] = inei+1; // >0 pour la réception
  }

  // On crée un accesseur pour la variable multi-env
  // Cette création est un peu couteuse (allocation)
  MultiEnvVar<DataType> menv_var(var_menv, m_mesh_material_mng);

  // On enchaine sur le device : 
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst

  // On remplit les buffers sur le DEVICE
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {

    auto byte_buf_snd_d = buf_snd_d.byteBuf(inei); // le buffer d'envoi pour inei sur le DEVICE
    auto byte_buf_snd_h = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei sur l'HOTE

    // On lit les valeurs de var pour les recopier dans le buffer d'envoi
    // byte_buf_snd_d = buffer dans lequel on va écrire
    // "buf_snd[inei] <= var_menv"
    async_pack_varmenv2buf(owned_evi_pn[inei], menv_var, byte_buf_snd_d, *(ref_queue.get()));

    // On enregistre un événement pour la fin de packing pour le voisin inei
    ref_queue->recordEvent(m_pack_events[inei]);

    // le transfert sur m_ref_queue_data ne pourra pas commencer
    // tant que l'événement m_pack_events[inei] ne sera pas arrivé
    m_ref_queue_data->waitEvent(m_pack_events[inei]);

    // transfert buf_snd_d[inei] => buf_snd_h[inei]
    async_transfer(byte_buf_snd_h, byte_buf_snd_d, *(m_ref_queue_data.get()));

    // On enregistre un événement pour la fin du transfert pour le voisin inei
    m_ref_queue_data->recordEvent(m_transfer_events[inei]);
  }
  
  // On amorce les envois sur l'HOTE dès qu'un transfert est terminé
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // Attente de la fin du transfert
    m_transfer_events[inei]->wait();

    // On amorce l'envoi
    auto byte_buf_snd = buf_snd_h.byteBuf(inei); // le buffer d'envoi pour inei
    requests[m_nb_nei+inei] = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    msg_types[m_nb_nei+inei] = -inei-1; // <0 pour l'envoi
  }

  // dans les faits, le stream est déjà synchronisé mais il faut terminer la RunQueue Arcane
  ref_queue->barrier();
  m_ref_queue_data->barrier();

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
        // on transfère les donnees de façon asynchrone pour le inei-ième voisin
        {
          auto byte_buf_rcv_h = buf_rcv_h.byteBuf(inei); // buffer des données reçues sur l'HOTE
          auto byte_buf_rcv_d = buf_rcv_d.byteBuf(inei); // buffer des données reçues à transférer sur le DEVICE

          // transfert buf_rcv_h[inei] => buf_rcv_d[inei]
          async_transfer(byte_buf_rcv_d, byte_buf_rcv_h, *(m_ref_queue_data.get()));

          // On enregistre un événement pour repérer la fin du transfert
          m_ref_queue_data->recordEvent(m_transfer_events[inei]);

          // Les kernels suivant sur ref_queue vont attendre l'occurence 
          // de l'événement m_transfer_events[inei] sur la queue m_ref_queue_data
          ref_queue->waitEvent(m_transfer_events[inei]);

          // Maintenant que byte_buf_rcv_d est sur DEVICE on peut enclencher le unpacking des données
          // "var_menv <= buf_rcv_d[inei]"
          async_unpack_buf2varmenv(ghost_evi_pn[inei], byte_buf_rcv_d, menv_var, *(ref_queue.get()));
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
  m_ref_queue_data->barrier();
  ref_queue->barrier();
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/
#define INST_VAR_SYNC_MNG_MULTI_MAT_SYNCHRONIZE(__DataType__) \
  template void VarSyncMng::multiMatSynchronize(Ref<RunQueue> ref_queue, CellMaterialVariableScalarRef<__DataType__> var)

INST_VAR_SYNC_MNG_MULTI_MAT_SYNCHRONIZE(Real);

