#ifndef MSG_PASS_COMPUTE_MAT_AND_SYNC_H
#define MSG_PASS_COMPUTE_MAT_AND_SYNC_H

#include <arcane/utils/ITraceMng.h>
#include <arcane/materials/MeshMaterialVariableSynchronizerList.h>

#include "msgpass/Algo1SyncDataMMatDH.h"
#include "msgpass/Algo1SyncDataMMatD.h"
#include "msgpass/VarSyncAlgo1.h"


/*---------------------------------------------------------------------------*/
/* Equivalent à un var.synchronize() où var est une variable multi-env       */ 
/*---------------------------------------------------------------------------*/
template<typename DataType>
void VarSyncMng::multiMatSynchronize(
    CellMaterialVariableScalarRef<DataType> var_menv, Ref<RunQueue> ref_queue, 
    eVarSyncVersion vs_version)
{
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var_menv);

  this->multiMatSynchronize(mvsl, ref_queue, vs_version);
}

/*---------------------------------------------------------------------------*/
/* Func func : traitement à appliquer sur un groupe de Cell                  */
/*                                                                           */
/* CellMaterialVariableScalarRef<DataType> var : variable multi-mat qui doit */
/* être synchronisée après les  calculs des items de bords                   */
/*                                                                           */
/* Applique le traitement func qui calcule var sur les items "own" de type   */
/* Cell et synchronise les items fantômes de var                             */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename Func, typename DataType>
void VarSyncMng::
computeMatAndSync(Func func, CellMaterialVariableScalarRef<DataType> var, eVarSyncVersion vs_version) {

  // Liste avec une seule variable
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var);

  computeMatAndSyncOnEvents<Func>(UniqueArray<Ref<ax::RunQueueEvent>>() /*liste vide d'événements*/,
      func, mvsl, vs_version);
}

/*---------------------------------------------------------------------------*/
/* depends_on_evts : liste d'événements dont va dépendre le traitement func  */
/*                                                                           */
/* Func func : traitement à appliquer sur un groupe de Cell                  */
/*                                                                           */
/* CellMaterialVariableScalarRef<DataType> var : variable multi-mat qui doit */
/* être synchronisée après les  calculs des items de bords                   */
/*                                                                           */
/* Applique le traitement func qui calcule var sur les items "own" de type   */
/* Cell et synchronise les items fantômes de var                             */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename Func, typename DataType>
void VarSyncMng::
computeMatAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, CellMaterialVariableScalarRef<DataType> var, eVarSyncVersion vs_version) {

  // Liste avec une seule variable
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var);

  computeMatAndSyncOnEvents<Func>(depends_on_evts,
      func, mvsl, vs_version);
}

/*---------------------------------------------------------------------------*/
/* Gestion avec une liste de variables multi-mat                             */
/*---------------------------------------------------------------------------*/
template<typename Func>
void VarSyncMng::
computeMatAndSync(Func func, MeshVariableSynchronizerList& vars, eVarSyncVersion vs_version) {
  computeMatAndSyncOnEvents<Func>(UniqueArray<Ref<ax::RunQueueEvent>>() /*liste vide d'événements*/,
      func, vars, vs_version);
}

/*---------------------------------------------------------------------------*/
/* depends_on_evts : liste d'événements dont va dépendre le traitement func  */
/*                                                                           */
/* Func func : traitement à appliquer sur un groupe de Cell                  */
/*                                                                           */
/* MeshVariableSynchronizerList& vars : liste des variables multi-mat qui    */
/*  doivent être synchronisées après les calculs des items de bords          */
/*                                                                           */
/* Applique le traitement func qui calcule les variables vars sur les items  */
/* "own" de type Cell et synchronise les items fantômes des vars             */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename Func>
void VarSyncMng::
computeMatAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, MeshVariableSynchronizerList& vars, eVarSyncVersion vs_version) {
  PROF_ACC_BEGIN(__FUNCTION__);

  using ItemType = Cell;

  // Pour qu'une queue attende des événements avant de commencer
  auto depends_on = [](Ref<RunQueue> rqueue, ArrayView<Ref<ax::RunQueueEvent>>& levts) {
    for(Integer iev=0 ; iev<levts.size() ; ++iev) {
      rqueue->waitEvent(levts[iev]);
    }
  };
  
  ITraceMng* tm = m_mesh->traceMng();
  if (vs_version == VS_bulksync_std ||
      vs_version == VS_bulksync_queue ||
      vs_version == VS_bulksync_evqueue) 
  {
    // Calcul sur tous les items "own"
    auto own_items = get_own_items<ItemType>(m_mesh);
    auto ref_queue = AcceleratorUtils::refQueueAsync(m_runner, QP_default); 

    // La queue va dépendre des événements de depends_on_evts
    depends_on(ref_queue, depends_on_evts);

    // Le calcul sur les items "propres" (hors items fantômes)
    func(own_items, ref_queue.get());
    ref_queue->barrier(); // on attend la fin du calcul sur tous les items "owns"

    // Puis comms
    PROF_ACC_BEGIN("syncVar");
    if (vs_version == VS_bulksync_std) {
      tm->debug() << "bulksync_std";
      MeshMaterialVariableSynchronizerList mmvsl(m_mesh_material_mng);
      for(auto v : vars.varsList()) {
        auto matv = v->materialVariable();
        if (matv) {
          mmvsl.add(matv);
        }
      }
      mmvsl.apply(); // les synchros regroupées en une
    } else if (vs_version == VS_bulksync_queue) {
      tm->debug() << "bulksync_queue";
      throw NotSupportedException(A_FUNCINFO,
        String::format("Invalid eVarSyncVersion={0}",(int)vs_version));
    } else {
      ARCANE_ASSERT(vs_version==VS_bulksync_evqueue,
          ("Ici, option differente de bulksync_evqueue"));
      tm->debug() << "bulksync_evqueue";
      Algo1SyncDataMMatDH sync_data_dh(vars, ref_queue, *m_a1_mmat_dh_pi);
      m_vsync_algo1->synchronize(&sync_data_dh);
    }
    PROF_ACC_END;

  } 
  else if (vs_version == VS_overlap_evqueue ||
      vs_version == VS_overlap_evqueue_d) 
  {
    SyncItems<ItemType>* sync_items = this->getSyncItems<ItemType>();

    // On amorce sur le DEVICE le calcul sur les items "own" sur le bord du
    // sous-domaine sur la queue m_ref_queue_bnd (_bnd = boundary)
    depends_on(m_ref_queue_bnd, depends_on_evts);
    func(sync_items->sharedItems(), m_ref_queue_bnd.get());
    // ici, le calcul n'est pas terminé sur le DEVICE

    // On amorce sur le DEVICE le calcul des items intérieurs dont ne dépendent
    // pas les comms sur la queue ref_queue_inr (_inr = inner)
    Ref<RunQueue> ref_queue_inr = AcceleratorUtils::refQueueAsync(m_runner, QP_default);
    depends_on(ref_queue_inr, depends_on_evts);
    func(sync_items->privateItems(), ref_queue_inr.get());

    // Sur la même queue de bord m_ref_queue_bnd, on amorce le packing des données
    // puis les comms MPI sur CPU, puis unpacking des données et on synchronise 
    // la queue m_ref_queue_bnd
    IAlgo1SyncData* sync_data=nullptr;
    if (vs_version == VS_overlap_evqueue) {
      tm->debug() << "overlap_evqueue";
      sync_data = new Algo1SyncDataMMatDH(vars, m_ref_queue_bnd, *m_a1_mmat_dh_pi);
    } else if (vs_version == VS_overlap_evqueue_d) {
      tm->debug() << "overlap_evqueue_d";
      sync_data = new Algo1SyncDataMMatD(vars, m_ref_queue_bnd, *m_a1_mmat_d_pi);
    }
    m_vsync_algo1->synchronize(sync_data);
    // ici, après cet appel, m_ref_queue_bnd est synchronisée
    delete sync_data;

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
    // Calcul sur tous les items "all"
    auto all_items = get_all_items<ItemType>(m_mesh);
    auto ref_queue = AcceleratorUtils::refQueueAsync(m_runner, QP_default); 
    depends_on(ref_queue, depends_on_evts);
    func(all_items, ref_queue.get());
    ref_queue->barrier();
  }
  PROF_ACC_END;
}

#endif

