#include "msgpass/VarSyncAlgo1.h"

#include <arcane/IParallelMng.h>
#include <arcane/utils/UniqueArray.h>
#include <arccore/base/FatalErrorException.h>

/*---------------------------------------------------------------------------*/
/* \class VarSyncAlgo1                                                       */
/* \brief Algorithm (v1) to synchronize mesh variables                       */
/*---------------------------------------------------------------------------*/

VarSyncAlgo1::VarSyncAlgo1(IParallelMng* pm, Int32ConstArrayView neigh_ranks) :
  m_pm          (pm),
  m_neigh_ranks (neigh_ranks)
{
  m_nb_nei = m_neigh_ranks.size();
}

VarSyncAlgo1::~VarSyncAlgo1() {
}

/*---------------------------------------------------------------------------*/
/* Synchronize variables encapsulated into sync_data                         */
/*---------------------------------------------------------------------------*/
void VarSyncAlgo1::synchronize(IAlgo1SyncData* sync_data)
{
  if (m_nb_nei==0 || sync_data->isEmpty()) {
    return;
  }

  // Step before the first communications
  sync_data->initComm();

  // L'échange proprement dit des valeurs de var
  UniqueArray<Parallel::Request> requests(2*m_nb_nei);
  IntegerUniqueArray msg_types(2*m_nb_nei); // nature des messages 

  // On amorce les réceptions
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    // On amorce la réception sur l'HOTE
    auto byte_buf_rcv = sync_data->recvBuf(inei); // le buffer de réception pour inei
    requests[inei] = m_pm->recv(byte_buf_rcv, rank_nei, /*blocking=*/false);
    msg_types[inei] = inei+1; // >0 pour la réception
  }

  sync_data->initSendings();
  
  // On amorce les envois sur l'HOTE dès qu'un transfert est terminé
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    Int32 rank_nei = m_neigh_ranks[inei]; // le rang du inei-ième voisin

    sync_data->finalizePackBeforeSend(inei);

    // On amorce l'envoi
    auto byte_buf_snd = sync_data->sendBuf(inei); // le buffer d'envoi pour inei
    requests[m_nb_nei+inei] = m_pm->send(byte_buf_snd, rank_nei, /*blocking=*/false);
    msg_types[m_nb_nei+inei] = -inei-1; // <0 pour l'envoi
  }

  sync_data->finalizeSendings();

  // Maitenant que toutes les requêtes de comms sont amorcées, il faut les terminer
  if (2*m_nb_nei!=requests.size()) {
    throw FatalErrorException(A_FUNCINFO, "Le nb de requetes n'est pas egal à 2 fois le nb de voisins");
  }

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

        nb_pending_rcv--; // on a une requete de reception en moins

        // On récupère l'indice du voisin
        Integer inei = pending_types[idone_req]-1;
        if (!(inei>=0 && inei<m_nb_nei)) {
          throw FatalErrorException(A_FUNCINFO, "Mauvais indice de voisin");
        }

        // Maintenant qu'on a reçu le buffer pour le inei-ième voisin, 
        // on post-traite les données reçues
        sync_data->unpackAfterRecv(inei);
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
    if (!(pending_requests.size()<=m_nb_nei)) {
      throw FatalErrorException(A_FUNCINFO,
          "Il ne peut pas rester un nb de requetes d'envoi supérieur au nb de voisins");
    }
    for(Integer msg_type : pending_types) {
      if (!(msg_type<0)) {
        throw FatalErrorException(A_FUNCINFO, 
            "Un message d'envoi doit avoir un type négatif ce qui n'est pas le cas");
      }
    }
    m_pm->waitAllRequests(pending_requests);
  }

  sync_data->finalizeReceipts();
}
