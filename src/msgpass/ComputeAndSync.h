#ifndef MSG_PASS_COMPUTE_AND_SYNC_H
#define MSG_PASS_COMPUTE_AND_SYNC_H

#include <arcane/utils/ITraceMng.h>

/*---------------------------------------------------------------------------*/
/* Soit ItemType =  MeshVariableRefT::ItemType                               */
/*                                                                           */
/* Func func : traitement à appliquer sur un groupe d'items ItemType         */
/*                                                                           */
/* MeshVariableRefT var : variable qui doit être synchronisée après les      */
/* calculs des items de bords                                                */
/*                                                                           */
/* Applique le traitement func qui calcule var sur les items "own" de type   */
/* ItemType et synchronise les items fantômes de var                         */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename Func, typename MeshVariableRefT>
void VarSyncMng::
computeAndSync(Func func, MeshVariableRefT var, eVarSyncVersion vs_version) {
  PROF_ACC_BEGIN(__FUNCTION__);

  using ItemType = typename MeshVariableRefT::ItemType;
  
  ITraceMng* tm = m_mesh->traceMng();
  if (vs_version == VS_bulksync_std ||
      vs_version == VS_bulksync_queue ||
      vs_version == VS_bulksync_evqueue) 
  {
    // Calcul sur tous les items "own"
    auto own_items = get_own_items<ItemType>(m_mesh);
    auto ref_queue = AcceleratorUtils::refQueueAsync(m_runner, QP_default); 
    // Le calcul sur les items "propres" (hors items fantômes)
    func(own_items, ref_queue.get());
    ref_queue->barrier(); // on attend la fin du calcul sur tous les items "owns"

    // Puis comms
    PROF_ACC_BEGIN("syncVar");
    if (vs_version == VS_bulksync_std) {
      tm->debug() << "bulksync_std";
      var.synchronize();
    } else if (vs_version == VS_bulksync_queue) {
      tm->debug() << "bulksync_queue";
      //this->globalSynchronize(var);
      this->globalSynchronizeQueue(ref_queue, var);
      //this->globalSynchronizeDevThr(var);
      //this->globalSynchronizeDevQueues(var);
    } else {
      ARCANE_ASSERT(vs_version==VS_bulksync_evqueue,
          ("Ici, option differente de bulksync_evqueue"));
      tm->debug() << "bulksync_evqueue";
      this->globalSynchronizeQueueEvent(ref_queue, var);
    }
    PROF_ACC_END;

  } 
  else if (vs_version == VS_overlap_evqueue ||
      vs_version == VS_overlap_evqueue_d) 
  {
    SyncItems<ItemType>* sync_items = this->getSyncItems<ItemType>();

    // On amorce sur le DEVICE le calcul sur les items "own" sur le bord du
    // sous-domaine sur la queue m_ref_queue_bnd (_bnd = boundary)
    func(sync_items->sharedItems(), m_ref_queue_bnd.get());
    // ici, le calcul n'est pas terminé sur le DEVICE

    // On amorce sur le DEVICE le calcul des items intérieurs dont ne dépendent
    // pas les comms sur la queue ref_queue_inr (_inr = inner)
    Ref<RunQueue> ref_queue_inr = AcceleratorUtils::refQueueAsync(m_runner, QP_default);
    func(sync_items->privateItems(), ref_queue_inr.get());

    // Sur la même queue de bord m_ref_queue_bnd, on amorce le packing des données
    // puis les comms MPI sur CPU, puis unpacking des données et on synchronise 
    // la queue m_ref_queue_bnd
    if (vs_version == VS_overlap_evqueue) {
      tm->debug() << "overlap_evqueue";
      this->globalSynchronizeQueueEvent(m_ref_queue_bnd, var);
    } else if (vs_version == VS_overlap_evqueue_d) {
      tm->debug() << "overlap_evqueue_d";
      this->globalSynchronizeQueueEventD(m_ref_queue_bnd, var);
    }
    // ici, après cet appel, m_ref_queue_bnd est synchronisée

    // On attend la terminaison des calculs intérieurs
    ref_queue_inr->barrier();
  } 
  else if (vs_version == VS_overlap_iqueue) 
  {
    tm->debug() << "overlap_iqueue";
    SyncItems<ItemType>* sync_items = this->getSyncItems<ItemType>();

    // On amorce sur le DEVICE le calcul sur les items "own" sur le bord du
    // sous-domaine sur la queue m_ref_queue_bnd (_bnd = boundary)
    func(sync_items->sharedItems(), m_ref_queue_bnd.get());
    // ici, le calcul n'est pas terminé sur le DEVICE
    
    // Sur la même queue de bord m_ref_queue_bnd, on amorce le packing des données
    auto ref_sync_req = this->iGlobalSynchronizeQueue(m_ref_queue_bnd, var);

    // On amorce sur le DEVICE le calcul des items intérieurs dont ne dépendent
    // pas les comms sur la queue ref_queue_inr (_inr = inner)
    Ref<RunQueue> ref_queue_inr = AcceleratorUtils::refQueueAsync(m_runner, QP_default);
    func(sync_items->privateItems(), ref_queue_inr.get());

    // Une fois le packing des données terminé sur m_ref_queue_bnd
    // on effectue les comms MPI sur CPU, 
    // puis unpacking des données et on synchronise la queue m_ref_queue_bnd
    ref_sync_req->wait();
    // ici, après cet appel, m_ref_queue_bnd est synchronisée

    // On attend la terminaison des calculs intérieurs
    ref_queue_inr->barrier();
  } 
  else
  {
    ARCANE_ASSERT(vs_version==VS_nosync,
        ("Ici, pas de synchro"));
    tm->debug() << "nosync";
    // Pas de synchro
    // Calcul sur tous les items "all"
    auto all_items = get_own_items<ItemType>(m_mesh);
    auto ref_queue = AcceleratorUtils::refQueueAsync(m_runner, QP_default); 
    func(all_items, ref_queue.get());
    ref_queue->barrier();
  }
  PROF_ACC_END;
}

#endif

