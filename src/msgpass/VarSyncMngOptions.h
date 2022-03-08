#ifndef MSG_PASS_VAR_SYNC_MNG_OPTIONS_H
#define MSG_PASS_VAR_SYNC_MNG_OPTIONS_H

/*! \brief Définit les implémentations de synchronisations d'une variable
 */
enum eVarSyncVersion {
  VS_nosync = 0,  // Pas de synchronisation
  VS_bulksync_std, // "Bulk-Synchronous" avec .synchronize() "classique" Arcane
  VS_bulksync_queue,  // "Bulk-Synchronous" avec packing/unpacking buf comm sur GPU
  VS_bulksync_evqueue, // "Bulk-Synchronous" avec packing/unpacking buf comm sur GPU avec events et waitSomeRequests
  VS_overlap_evqueue, // Recouvrement items shared+packing/unpacking GPU+comms (en utilisant des events) par calculs itemss private
  VS_overlap_evqueue_d, // Idem que VS_overlap_evqueue mais comms avec adresses DEVICE (GPU-aware)
  VS_overlap_iqueue // Recouvrement : traitements shared et private concurrents et asynchrones + iGlobalSynchronizeQueue
};

#endif

