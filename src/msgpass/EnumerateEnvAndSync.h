#ifndef MSG_PASS_ENUMERATE_ENV_AND_SYNC_H
#define MSG_PASS_ENUMERATE_ENV_AND_SYNC_H

#include <arcane/utils/ITraceMng.h>
#include <arcane/materials/MeshMaterialVariableSynchronizerList.h>

#include "msgpass/VarSyncAlgo1.h"


/*---------------------------------------------------------------------------*/
/* Func func : traitement à appliquer sur un groupe de Cell                  */
/*                                                                           */
/* CellMaterialVariableScalarRef<DataType> var : variable multi-mat qui doit */
/* être synchronisée après les calculs des EnvVarIndex(es) de bords          */
/*                                                                           */
/* Applique le traitement func qui calcule var sur les EnvVarIndex(es) "own" */
/* et synchronise les EnvVarIndex(es) fantômes de var                        */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename Func, typename DataType>
void VarSyncMng::
enumerateEnvAndSync(Func func, CellMaterialVariableScalarRef<DataType> var, eVarSyncVersion vs_version) {

  // Liste avec une seule variable
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var);

  enumerateEnvAndSyncOnEvents<Func>(UniqueArray<Ref<ax::RunQueueEvent>>() /*liste vide d'événements*/,
      func, mvsl, vs_version);
}

/*---------------------------------------------------------------------------*/
/* depends_on_evts : liste d'événements dont va dépendre le traitement func  */
/*                                                                           */
/* Func func : traitement à appliquer sur un groupe de Cell                  */
/*                                                                           */
/* CellMaterialVariableScalarRef<DataType> var : variable multi-mat qui doit */
/* être synchronisée après les calculs des EnvVarIndex(es) de bords          */
/*                                                                           */
/* Applique le traitement func qui calcule var sur les EnvVarIndex(es) "own" */
/* et synchronise les EnvVarIndex(es) fantômes de var                        */
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename Func, typename DataType>
void VarSyncMng::
enumerateEnvAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, CellMaterialVariableScalarRef<DataType> var, eVarSyncVersion vs_version) {

  // Liste avec une seule variable
  MeshVariableSynchronizerList mvsl(this);
  mvsl.add(var);

  enumerateEnvAndSyncOnEvents<Func>(depends_on_evts,
      func, mvsl, vs_version);
}

/*---------------------------------------------------------------------------*/
/* Gestion avec une liste de variables multi-mat                             */
/*---------------------------------------------------------------------------*/
template<typename Func>
void VarSyncMng::
enumerateEnvAndSync(Func func, MeshVariableSynchronizerList& vars, eVarSyncVersion vs_version) {
  enumerateEnvAndSyncOnEvents<Func>(UniqueArray<Ref<ax::RunQueueEvent>>() /*liste vide d'événements*/,
      func, vars, vs_version);
}

/*---------------------------------------------------------------------------*/
/* depends_on_evts : liste d'événements dont va dépendre le traitement func  */
/*                                                                           */
/* Func func : traitement à appliquer sur chaque environnement               */
/*                                                                           */
/* MeshVariableSynchronizerList& vars : liste des variables multi-mat qui    */
/* doivent être synchronisées après les calculs des EnvVarIndex(es) de bords */
/*                                                                           */
/* Applique le traitement func qui calcule les variables vars sur les        */
/* EnvVarIndex(es) "own" et synchronise les EnvVarIndex(es) fantômes des vars*/
/* Note : ceci peut se faire dans n'importe quel ordre (vs_version)          */
/*---------------------------------------------------------------------------*/
template<typename Func>
void VarSyncMng::
enumerateEnvAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, MeshVariableSynchronizerList& vars, eVarSyncVersion vs_version) {
  PROF_ACC_BEGIN(__FUNCTION__);

  // Pour qu'une queue attende des événements avant de commencer
  auto depends_on = [](Ref<RunQueue> rqueue, ArrayView<Ref<ax::RunQueueEvent>>& levts) {
    for(Integer iev=0 ; iev<levts.size() ; ++iev) {
      rqueue->waitEvent(levts[iev]);
    }
  };

  // Pour que l'hôte attende des événements avant de commencer
  auto wait_for = [](ArrayView<Ref<ax::RunQueueEvent>>& levts) {
    for(Integer iev=0 ; iev<levts.size() ; ++iev) {
      levts[iev]->wait();
    }
  };

  // Pour appliquer func pour tous les environnements
  auto allenv_func = [&](ConstMultiArray2View<EnvVarIndex> levis_penv,
      MultiAsyncRunQueue* menv_queue)
  {
    // chaque environnement est calculé en parallèle
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();

      // calcul asynchrone car menv_queue->queue(env_id) est asynchrone
      func(env, levis_penv[env_id], &(menv_queue->queue(env_id)));
      // TODO : enregistrer un événement sur menv_queue->queue(env_id)
    }
  };
  
  ITraceMng* tm = m_mesh->traceMng();
  if (vs_version == VS_bulksync_std ||
      vs_version == VS_bulksync_queue ||
      vs_version == VS_bulksync_evqueue ||
      vs_version == VS_bulksync_evqueue_d) 
  {
    // Calcul sur les EnvVarIndex(es) intérieurs "own"
    auto own_levis_penv = m_sync_evi->ownEviPenv();

    // La queue va dépendre des événements de depends_on_evts
    wait_for(depends_on_evts);

    // Le calcul sur les EnvVarIndex(es) "propres" (hors items fantômes)
    allenv_func(own_levis_penv, m_menv_queue_inr);

    // on attend la fin du calcul sur tous les environnements
    m_menv_queue_inr->waitAllQueues();

    // Puis comms
    // On utilise m_ref_queue_inr qui a une priorité par défault
    this->synchronize(vars, m_ref_queue_inr, vs_version);

  } 
  else if (vs_version == VS_overlap_evqueue ||
      vs_version == VS_overlap_evqueue_d) 
  {
    // On amorce sur le DEVICE le calcul sur les EnvVarIndex(es) "shared" sur 
    // le bord du sous-domaine sur la queue m_menv_queue_bnd (_bnd = boundary)
    auto shared_levis_penv = m_sync_evi->sharedEviPenv();

    // La queue va dépendre des événements de depends_on_evts
    wait_for(depends_on_evts);

    // Le calcul sur les EnvVarIndex(es) "shared" (bord intérieur du sous-domaine)
    // sur des queues prioritaires m_menv_queue_bnd
    allenv_func(shared_levis_penv, m_menv_queue_bnd);

    // On amorce sur le DEVICE le calcul des EnvVarIndex(es) intérieurs dont ne dépendent
    // pas les comms sur la queue m_menv_queue_inr (_inr = inner)
    auto private_levis_penv = m_sync_evi->privateEviPenv();
    allenv_func(private_levis_penv, m_menv_queue_inr);

    // on attend la fin du calcul "shared" sur tous les environnements 
    // avant d'amorcer les comms
    m_menv_queue_bnd->waitAllQueues();
    // Sur la queue prioritaire m_ref_queue_bnd, on amorce le packing des données
    // puis les comms MPI sur CPU, puis unpacking des données et on synchronise 
    this->synchronize(vars, m_ref_queue_bnd, vs_version);

    // On attend la terminaison des calculs intérieurs
    m_menv_queue_inr->waitAllQueues();
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
    // Calcul sur tous les EnvVarIndex(es) "all"
    auto all_levis_penv = m_sync_evi->allEviPenv();

    // La queue va dépendre des événements de depends_on_evts
    wait_for(depends_on_evts);

    // Le calcul sur tous les EnvVarIndex(es) (fantômes compris)
    allenv_func(all_levis_penv, m_menv_queue_inr);

    // on attend la fin du calcul sur tous les environnements
    m_menv_queue_inr->waitAllQueues();
  }
  PROF_ACC_END;
}

#endif

