#include "msgpass/Algo1SyncDataGlobDH.h"
#include "msgpass/PackTransfer.h"

Algo1SyncDataGlobDH::PersistentInfo::PersistentInfo(
    Integer nb_nei,
    Runner& runner,
    SyncBuffers* sync_buffers) :
  m_sync_buffers (sync_buffers),
  m_nb_nei       (nb_nei)
{
  m_pack_events.resize(m_nb_nei);
  m_transfer_events.resize(m_nb_nei);
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    m_pack_events[inei]     = makeEventRef(runner);
    m_transfer_events[inei] = makeEventRef(runner);
  }

  m_ref_queue_data = AcceleratorUtils::refQueueAsync(runner, QP_high);
}

Algo1SyncDataGlobDH::PersistentInfo::~PersistentInfo()
{
}

/*---------------------------------------------------------------------------*/
/* \class Algo1SyncDataGlobDH                                                */
/* \brief Implementation of IAlgo1SyncData for global variables (ie non multi-mat) */
/*   packing/unpacking on Device                                             */
/*   communicating (MPI) on Host                                             */
/*---------------------------------------------------------------------------*/

Algo1SyncDataGlobDH::Algo1SyncDataGlobDH(
    MeshVariableSynchronizerList& vars,
    Ref<RunQueue> ref_queue,
    Algo1SyncDataGlobDH::PersistentInfo& pi
    ) :
  m_vars         (vars),
  m_ref_queue    (ref_queue),
  m_pi           (pi)
{
}

Algo1SyncDataGlobDH::~Algo1SyncDataGlobDH() {
}

/*---------------------------------------------------------------------------*/
/* True if there is no data to synchronize                                   */
/*---------------------------------------------------------------------------*/
bool Algo1SyncDataGlobDH::isEmpty() {
  return m_vars.varsList().size()==0;
}

/*---------------------------------------------------------------------------*/
/* Initialization step before beginning communications                       */
/*---------------------------------------------------------------------------*/
void Algo1SyncDataGlobDH::initComm() {

  auto lvars = m_vars.varsList();

  // Asynchronous pointers tranfer onto device
  m_vars.asyncHToD(*(m_ref_queue.get()));

  // On prévoit une taille max du buffer qui va contenir tous les messages
  Int64 buf_estim_sz=0;
  for(auto var : lvars) {
    buf_estim_sz += var->estimatedMaxBufSz();
  }

  m_pi.m_sync_buffers->resetBuf();
  // Le buffer de tous les messages est réalloué si pas assez de place
  m_pi.m_sync_buffers->allocIfNeeded(buf_estim_sz);

  Integer nb_nei = m_pi.m_nb_nei;
  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur l'HOTE (_h et "0")
  m_buf_snd_h = m_pi.m_sync_buffers->multiBufViewVars(lvars, nb_nei, IMeshVarSync::IS_owned, 0);
  m_buf_rcv_h = m_pi.m_sync_buffers->multiBufViewVars(lvars, nb_nei, IMeshVarSync::IS_ghost, 0);

  // On récupère les adresses et tailles des buffers d'envoi et de réception 
  // sur le DEVICE (_d et "1")
  m_buf_snd_d = m_pi.m_sync_buffers->multiBufViewVars(lvars, nb_nei, IMeshVarSync::IS_owned, 1);
  m_buf_rcv_d = m_pi.m_sync_buffers->multiBufViewVars(lvars, nb_nei, IMeshVarSync::IS_ghost, 1);
}


/*---------------------------------------------------------------------------*/
/* Get the receive buffer for the neighbour inei                             */
/*---------------------------------------------------------------------------*/
ArrayView<Byte> Algo1SyncDataGlobDH::recvBuf(Integer inei) {
  return m_buf_rcv_h.multiView(inei).rangeView();
}

/*---------------------------------------------------------------------------*/
/* Get the send buffer for the neighbour inei                                */
/*---------------------------------------------------------------------------*/
ArrayView<Byte> Algo1SyncDataGlobDH::sendBuf(Integer inei) {
  return m_buf_snd_h.multiView(inei).rangeView();
}

/*---------------------------------------------------------------------------*/
/* Prepare the sendings for every neighbour                                  */
/*---------------------------------------------------------------------------*/
void Algo1SyncDataGlobDH::initSendings() {

  auto lvars = m_vars.varsList();
  Integer nb_var = lvars.size();

  // On enchaine sur le device : 
  //    copie de var_dev dans buf_dev 
  //    puis transfert buf_dev => buf_hst

  // On remplit les buffers sur le DEVICE
  for(Integer inei=0 ; inei<m_pi.m_nb_nei ; ++inei) {

    auto byte_buf_snd_d = m_buf_snd_d.multiView(inei); // le buffer d'envoi pour inei sur le DEVICE
    auto byte_buf_snd_h = m_buf_snd_h.multiView(inei); // le buffer d'envoi pour inei sur l'HOTE

    for(Integer ivar=0 ; ivar<nb_var ; ++ivar) {
      // On lit les valeurs de lvars[ivar] pour les recopier dans le buffer d'envoi
      // byte_buf_var_d = buffer dans lequel on va écrire
      // "buf_snd[inei] <= var_menv"
      auto byte_buf_var_d = byte_buf_snd_d.byteBuf(ivar);
      lvars[ivar]->asyncPackOwnedIntoBuf(inei, byte_buf_var_d, *(m_ref_queue.get()));
    }

    // On enregistre un événement pour la fin de packing pour le voisin inei
    m_ref_queue->recordEvent(m_pi.m_pack_events[inei]);

    // le transfert sur m_pi.m_ref_queue_data ne pourra pas commencer
    // tant que l'événement m_pi.m_pack_events[inei] ne sera pas arrivé
    m_pi.m_ref_queue_data->waitEvent(m_pi.m_pack_events[inei]);

    // transfert m_buf_snd_d[inei] => m_buf_snd_h[inei]
    async_transfer(byte_buf_snd_h, byte_buf_snd_d, *(m_pi.m_ref_queue_data.get()));

    // On enregistre un événement pour la fin du transfert pour le voisin inei
    m_pi.m_ref_queue_data->recordEvent(m_pi.m_transfer_events[inei]);
  }
}

/*---------------------------------------------------------------------------*/
/* Prepare the sending for one neighbour                                     */
/*---------------------------------------------------------------------------*/
void Algo1SyncDataGlobDH::finalizePackBeforeSend(Integer inei) {
  // Attente de la fin du transfert
  m_pi.m_transfer_events[inei]->wait();
}

/*---------------------------------------------------------------------------*/
/* Finalize the sendings for every neighbour                                 */
/*---------------------------------------------------------------------------*/
void Algo1SyncDataGlobDH::finalizeSendings() {
  // dans les faits, le stream est déjà synchronisé mais il faut terminer la RunQueue Arcane
  m_ref_queue->barrier();
  m_pi.m_ref_queue_data->barrier();
}

/*---------------------------------------------------------------------------*/
/* Data treatement after receiving the message from the neighbour inei       */
/*---------------------------------------------------------------------------*/
void Algo1SyncDataGlobDH::unpackAfterRecv(Integer inei) {

  auto lvars = m_vars.varsList();
  Integer nb_var = lvars.size();

  // Maintenant qu'on a reçu le buffer pour le inei-ième voisin, 
  // on transfère les donnees de façon asynchrone pour le inei-ième voisin

  auto byte_buf_rcv_h = m_buf_rcv_h.multiView(inei); // buffer des données reçues sur l'HOTE
  auto byte_buf_rcv_d = m_buf_rcv_d.multiView(inei); // buffer des données reçues à transférer sur le DEVICE

  // transfert m_buf_rcv_h[inei] => m_buf_rcv_d[inei]
  async_transfer(byte_buf_rcv_d, byte_buf_rcv_h, *(m_pi.m_ref_queue_data.get()));

  // On enregistre un événement pour repérer la fin du transfert
  m_pi.m_ref_queue_data->recordEvent(m_pi.m_transfer_events[inei]);

  // Les kernels suivants sur m_ref_queue vont attendre l'occurence 
  // de l'événement m_pi.m_transfer_events[inei] sur la queue m_pi.m_ref_queue_data
  m_ref_queue->waitEvent(m_pi.m_transfer_events[inei]);

  // Maintenant que byte_buf_rcv_d est sur DEVICE on peut enclencher le unpacking des données
  // variable par variable
  // "var_menv <= m_buf_rcv_d[inei]"
  for(Integer ivar=0 ; ivar<nb_var ; ++ivar) {
    auto byte_buf_var_d = byte_buf_rcv_d.byteBuf(ivar);
    lvars[ivar]->asyncUnpackGhostFromBuf(inei, byte_buf_var_d, *(m_ref_queue.get()));
  }
}

/*---------------------------------------------------------------------------*/
/* Finalize the receipts for every neighbour                                 */
/*---------------------------------------------------------------------------*/
void Algo1SyncDataGlobDH::finalizeReceipts() {
  // on attend la terminaison de tous les unpacks asynchrones
  m_pi.m_ref_queue_data->barrier();
  m_ref_queue->barrier();
}

/*---------------------------------------------------------------------------*/
/* Finalize if no communication needed                                       */
/*---------------------------------------------------------------------------*/
void Algo1SyncDataGlobDH::finalizeWoComm() {
  m_ref_queue->barrier();
}

