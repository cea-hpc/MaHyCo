#ifndef MSG_PASS_SYNC_AND_COMPUTE_H
#define MSG_PASS_SYNC_AND_COMPUTE_H

#include <arcane/utils/ITraceMng.h>

/*---------------------------------------------------------------------------*/
/* Func func : traitement à appliquer sur un groupe de ItemType              */
/*                                                                           */
/* CellMaterialVariableScalarRef<DataType> var : variable multi-mat qui doit */
/* être synchronisée après les  calculs des items de bords                   */
/*                                                                           */
/* Synchronise les items fantômes de var puis applique le traitement func    */
/* qui calcule sur les items item_group de type ItemType                     */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename Func, typename DataType>
void VarSyncMng::
syncAndCompute(CellMaterialVariableScalarRef<DataType> var, 
    ItemGroupT<ItemType> item_group, Func func, 
    eVarSyncVersion vs_version) {

  // Liste avec une seule variable
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var);

  syncAndComputeOnEvents<ItemType, Func>(
      UniqueArray<Ref<ax::RunQueueEvent>>() /*liste vide d'événements*/,
      mvsl, item_group, func, vs_version);
}

/*---------------------------------------------------------------------------*/
/* depends_on_evts : liste d'événements dont va dépendre le traitement func  */
/*                                                                           */
/* Func func : traitement à appliquer sur un groupe de ItemType              */
/*                                                                           */
/* CellMaterialVariableScalarRef<DataType> var : variable multi-mat qui doit */
/* être synchronisée après les  calculs des items de bords                   */
/*                                                                           */
/* Synchronise les items fantômes de var puis applique le traitement func    */
/* qui calcule sur les items item_group de type ItemType                     */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename Func, typename DataType>
void VarSyncMng::
syncAndComputeOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    CellMaterialVariableScalarRef<DataType> var, 
    ItemGroupT<ItemType> item_group, Func func, 
    eVarSyncVersion vs_version) {

  // Liste avec une seule variable
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var);

  syncAndComputeOnEvents<item_group, Func>(depends_on_evts,
      mvsl, item_group, func, vs_version);
}

/*---------------------------------------------------------------------------*/
/*  TODO                                                                     */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename Func, typename MeshVariableRefT>
void VarSyncMng::
syncAndCompute(MeshVariableRefT var, 
    ItemGroupT<ItemType> item_group, Func func, 
    eVarSyncVersion vs_version) {

  // Liste avec une seule variable
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var);

  syncAndCompute<ItemType, Func>(mvsl, item_group, func, vs_version);
}

/*---------------------------------------------------------------------------*/
/* Gestion avec une liste de variables multi-mat                             */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename Func>
void VarSyncMng::
syncAndCompute(MeshVariableSynchronizerList& vars, 
    ItemGroupT<ItemType> item_group, Func func,
    eVarSyncVersion vs_version) {

  syncAndComputeOnEvents<ItemType, Func>(
      UniqueArray<Ref<ax::RunQueueEvent>>() /*liste vide d'événements*/,
      vars, item_group, func, vs_version);
}

/*---------------------------------------------------------------------------*/
/* depends_on_evts : liste d'événements dont va dépendre le traitement func  */
/*                                                                           */
/* Func func : traitement à appliquer sur un groupe de ItemType              */
/*                                                                           */
/* MeshVariableSynchronizerList& vars : liste des variables globales et/ou   */
/*  multi-mat qui doivent être synchronisées AVANT les calculs des items de  */
/*  bords                                                                    */
/*                                                                           */
/* Synchronise les items fantômes des vars puis applique le traitement func  */
/* qui calcule sur les items item_group de type ItemType                     */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename Func>
void VarSyncMng::
syncAndComputeOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    MeshVariableSynchronizerList& vars, 
    ItemGroupT<ItemType> item_group, Func func, 
    eVarSyncVersion vs_version) 
{
  PROF_ACC_BEGIN(__FUNCTION__);

  // Les seuls groupes supportés sont "own" ou "all"
  typename SyncItems<ItemType>::eGroupCategory group_category;
  if (item_group == get_own_items<ItemType>(m_mesh)) {
    group_category = SyncItems<ItemType>::GC_own;
  } else if (item_group == get_all_items<ItemType>(m_mesh)) {
    group_category = SyncItems<ItemType>::GC_all;
  } else {
    throw NotSupportedException(A_FUNCINFO, 
	String("Les seuls groupes supportés sont \"own\" ou \"all\""));
  }

  // Pour qu'une queue attende des événements avant de commencer
  auto depends_on = [](Ref<RunQueue> rqueue, ArrayView<Ref<ax::RunQueueEvent>>& levts) {
    for(Integer iev=0 ; iev<levts.size() ; ++iev) {
      rqueue->waitEvent(levts[iev]);
    }
  };

  if (vs_version == VS_auto) {
    vs_version = defaultVarSyncVersion();
  }
  
  ITraceMng* tm = m_mesh->traceMng();
  if (vs_version == VS_bulksync_std ||
      vs_version == VS_bulksync_queue ||
      vs_version == VS_bulksync_evqueue) 
  {
    auto ref_queue = AcceleratorUtils::refQueueAsync(m_runner, QP_default); 

    // La queue va dépendre des événements de depends_on_evts
    depends_on(ref_queue, depends_on_evts);

    // Comms
    this->synchronize(vars, ref_queue, vs_version);

    // Puis calcul sur les items demandés
    func(item_group, ref_queue.get());
    ref_queue->barrier(); // on attend la fin du calcul sur tous les items demandés
  } 
  else if (vs_version == VS_overlap_evqueue ||
      vs_version == VS_overlap_evqueue_d) 
  {
    SyncItems<ItemType>* sync_items = this->getSyncItems<ItemType>();

    // On amorce sur le DEVICE le calcul des items intérieurs dont ne dépendent
    // pas les comms sur la queue ref_queue_inr (_inr = inner)
    Ref<RunQueue> ref_queue_inr = AcceleratorUtils::refQueueAsync(m_runner, QP_default);
    depends_on(ref_queue_inr, depends_on_evts);
    func(sync_items->privateItems(), ref_queue_inr.get());
    // ici, le calcul sur ref_queue_inr n'est pas terminé sur le DEVICE

    // Sur la queue de bord m_ref_queue_bnd, on amorce le packing des données
    // puis les comms MPI, puis unpacking des données et on synchronise 
    // la queue m_ref_queue_bnd
    depends_on(m_ref_queue_bnd, depends_on_evts);
    this->synchronize(vars, m_ref_queue_bnd, vs_version);

    // Puis on amorce sur le DEVICE le calcul sur les items dont dépendent
    // les comms sur la queue m_ref_queue_bnd (_bnd = boundary)
    func(sync_items->boundaryItems(group_category), m_ref_queue_bnd.get());
    m_ref_queue_bnd->barrier();

    // On attend la terminaison des calculs intérieurs
    ref_queue_inr->barrier();
  } 
  else if (vs_version == VS_overlap_iqueue) 
  {
    tm->debug() << "overlap_iqueue";
    throw NotSupportedException(A_FUNCINFO,
        String::format("Invalid eVarSyncVersion={0}",(int)vs_version));
  } 
  else
  {
    ARCANE_ASSERT(vs_version==VS_nosync,
        ("Ici, pas de synchro"));
    tm->debug() << "nosync";
    // Pas de synchro
    // Calcul sur tous les items demandés
    auto ref_queue = AcceleratorUtils::refQueueAsync(m_runner, QP_default); 
    depends_on(ref_queue, depends_on_evts);
    func(item_group, ref_queue.get());
    ref_queue->barrier();
  }
  PROF_ACC_END;
}

#endif

