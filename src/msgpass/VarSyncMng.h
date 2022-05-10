#ifndef VAR_SYNC_MNG_H
#define VAR_SYNC_MNG_H

#include "accenv/AcceleratorUtils.h"
#include <arcane/IMesh.h>
#include <arcane/materials/IMeshMaterialMng.h>
#include <arcane/materials/MeshMaterialVariableRef.h>
#include <arcane/utils/MultiArray2.h>
#include <arcane/IParallelMng.h>
#include <arcane/accelerator/core/RunQueueEvent.h>

#include "msgpass/SyncItems.h"
#include "msgpass/SyncEnvIndexes.h"
#include "msgpass/SyncBuffers.h"
#include "msgpass/VarSyncMngOptions.h"
#include "msgpass/MeshVariableSynchronizerList.h"
#include "msgpass/VarSyncAlgo1.h"
#include "msgpass/Algo1SyncDataMMatDH.h"
#include "msgpass/Algo1SyncDataMMatD.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Encapsule une requête de synchronisation sur une variable globale         */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
class GlobalSyncRequest {
 public:
  GlobalSyncRequest(
    SyncItems<ItemType>* sync_items,
    SyncBuffers* sync_buffers,
    IParallelMng* pm,
    Int32ConstArrayView neigh_ranks,
    Ref<RunQueue> ref_queue,
    MeshVarRefT<ItemType, DataType> var);

  virtual ~GlobalSyncRequest();

  // Termine les comms de m_requests, synchronise la queue m_ref_queue, 
  // met à jour les items fantômes de var
  void wait();

 protected:
  bool m_is_over=false;  //! Indique si la requête est terminée
  IParallelMng* m_pm=nullptr;  //! pour effectuer les send/receive proprement dit
  Int32ConstArrayView m_neigh_ranks;  //! Liste des rangs des voisins
  UniqueArray<Parallel::Request> m_requests;  //! les requetes send/receive proprement dites

  ConstMultiArray2View<Integer> m_ghost_item_idx_pn;  //! par voisin, items fantômes à mettre à jour
  
  // Les buffers pour tous les voisins envois/réceptions sur host et device
  MultiBufView m_buf_snd_h;
  MultiBufView m_buf_snd_d;
  MultiBufView m_buf_rcv_h;
  MultiBufView m_buf_rcv_d;

  Ref<RunQueue> m_ref_queue;  //! la queue sur laquelle effectuer les opérations sur GPU
  MeshVarRefT<ItemType, DataType> m_var;  //! variable Arcane dont il faut mettre les items fantômes à jour
};

/*---------------------------------------------------------------------------*/
/* Gère les synchronisations des mailles fantômes par Message Passing        */
/*---------------------------------------------------------------------------*/
class VarSyncMng {
 public:
  VarSyncMng(IMesh* mesh, ax::Runner& runner, AccMemAdviser* acc_mem_adv);
  virtual ~VarSyncMng();

  //! Initialise les futures synchros multi-mat
  void initSyncMultiEnv(IMeshMaterialMng* mesh_material_mng);

  //! Remet à jour les synchros multi-mat quand la carte des environnements a changé
  void updateSyncMultiEnv();

  //! Retourne vrai si un GPU est dispo pour exécuter les calculs
  bool isAcceleratorAvailable() const;

  //! Retourne vrai si on peut utiliser les adresses dans DEVICE pour les comms
  bool isDeviceAware() const;

  //! Buffer d'adresses pour gérer les côuts des allocations
  BufAddrMng* bufAddrMng();

  //! Pour gérer les EnvVarIndex(es) pour les comms
  SyncEnvIndexes* syncEnvIndexes();

  // Retourne l'instance de SyncItems<T> en fonction de T
  template<typename ItemType>
  SyncItems<ItemType>* getSyncItems();

  // Equivalent à un var.synchronize() où var est une variable globale (i.e. non multi-mat)
  template<typename MeshVariableRefT>
  void globalSynchronize(MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  // La queue asynchrone ref_queue est synchronisé en fin d'appel
  template<typename MeshVariableRefT>
  void globalSynchronizeQueue(Ref<RunQueue> ref_queue, MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  template<typename MeshVariableRefT>
  void globalSynchronizeDevThr(MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  template<typename MeshVariableRefT>
  void globalSynchronizeDevQueues(MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice
  template<typename MeshVariableRefT>
  void globalSynchronizeQueueEvent(Ref<RunQueue> ref_queue, MeshVariableRefT var);

  // Equivalent à un globalSynchronize pour lequel les données de var sont sur le DEVice et les comms se font avec les adresses du DEVICE
  template<typename MeshVariableRefT>
  void globalSynchronizeQueueEventD(Ref<RunQueue> ref_queue, MeshVariableRefT var);

  // Amorce un globalSynchronizeQueue
  template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
  Ref<GlobalSyncRequest<ItemType, DataType, MeshVarRefT> > iGlobalSynchronizeQueue(Ref<RunQueue> ref_queue, MeshVarRefT<ItemType, DataType> var);

  // Overlapping entre calcul et communications pour variable "globale"
  template<typename Func, typename MeshVariableRefT>
  void computeAndSync(Func func, MeshVariableRefT var, eVarSyncVersion vs_version=VS_overlap_evqueue);


  //! Equivalent à un var.synchronize() où var est une variable multi-mat
  template<typename DataType>
  void multiMatSynchronize(CellMaterialVariableScalarRef<DataType> var_menv, Ref<RunQueue> ref_queue, 
      eVarSyncVersion vs_version=VS_bulksync_evqueue);

  //! Maj des mailles fantômes d'une liste de variables multi-mat
  void multiMatSynchronize(MeshVariableSynchronizerList& vars, Ref<RunQueue> ref_queue, 
      eVarSyncVersion vs_version=VS_bulksync_evqueue);

  /* computeMatAndSync[OnEvents] */

  // Overlapping entre calcul et communications pour variable multi-mat
  template<typename Func, typename DataType>
  void computeMatAndSync(Func func, CellMaterialVariableScalarRef<DataType> var, 
      eVarSyncVersion vs_version=VS_overlap_evqueue);

  //! Une fois déroulés les événements de depends_on_evts, 
  // overlapping entre calcul et communications pour variable multi-mat
  template<typename Func, typename DataType>
  void computeMatAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, CellMaterialVariableScalarRef<DataType> var, 
    eVarSyncVersion vs_version=VS_overlap_evqueue);

  // Overlapping entre calcul et communications pour une liste de variables multi-mat
  template<typename Func>
  void computeMatAndSync(Func func, MeshVariableSynchronizerList& vars, 
      eVarSyncVersion vs_version=VS_overlap_evqueue);

  //! Une fois déroulés les événements de depends_on_evts, 
  // overlapping entre calcul et communications pour une liste de variables multi-mat
  template<typename Func>
  void computeMatAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
      Func func, MeshVariableSynchronizerList& vars, 
      eVarSyncVersion vs_version=VS_overlap_evqueue);

  /* enumerateEnvAndSync[OnEvents] */

  template<typename Func, typename DataType>
  void enumerateEnvAndSync(Func func, CellMaterialVariableScalarRef<DataType> var, 
      eVarSyncVersion vs_version=VS_overlap_evqueue);

  template<typename Func, typename DataType>
  void enumerateEnvAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, CellMaterialVariableScalarRef<DataType> var, 
    eVarSyncVersion vs_version=VS_overlap_evqueue);

  template<typename Func>
  void enumerateEnvAndSync(Func func, MeshVariableSynchronizerList& vars, 
      eVarSyncVersion vs_version=VS_overlap_evqueue);

  template<typename Func>
  void enumerateEnvAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, MeshVariableSynchronizerList& vars, 
    eVarSyncVersion vs_version=VS_overlap_evqueue);

 protected:
  
  // Pré-allocation des buffers de communication pour miniser le nb de réallocations
  void _preAllocBuffers();

 protected:

  IMesh* m_mesh=nullptr;
  IMeshMaterialMng* m_mesh_material_mng=nullptr;
  AccMemAdviser* m_acc_mem_adv=nullptr;

  ax::Runner& m_runner;

  bool m_is_device_aware=false; //! Vrai si l'on peut effectuer les comms avec adresses sur GPU

  IParallelMng* m_pm;  //! pour effectuer les send/receive proprement dit

  Integer m_nb_nei;  //! Nb de voisins (m_neigh_ranks.size())
  Int32ConstArrayView m_neigh_ranks;  //! Liste des rangs des voisins

  SyncItems<Cell>* m_sync_cells=nullptr;
  SyncItems<Node>* m_sync_nodes=nullptr;

  SyncBuffers* m_sync_buffers=nullptr;

  MultiAsyncRunQueue* m_neigh_queues=nullptr;  //! Pour gérer plusieurs queues pour les voisins

  Ref<ax::RunQueue> m_ref_queue_inr;  //! Référence sur queue standard pour traitement items intérieurs
  Ref<ax::RunQueue> m_ref_queue_bnd;  //! Référence sur queue prioritaire pour traitement items bords
  Ref<ax::RunQueue> m_ref_queue_data;  //! Référence sur une queue prioritaire pour le transfert des données
  UniqueArray<Ref<ax::RunQueueEvent>> m_pack_events;  //! Les evenements pour le packing des données
  UniqueArray<Ref<ax::RunQueueEvent>> m_transfer_events;  //! Les evenements pour le transfert des données

  SyncEnvIndexes* m_sync_evi=nullptr;  //! Pour gérer les EnvVarIndex(es) pour les comms

  MultiAsyncRunQueue* m_menv_queue_inr=nullptr;  //!< queues pour traiter les environnements de façon asynchrone à l'intérieur du sous-domaine
  MultiAsyncRunQueue* m_menv_queue_bnd=nullptr;  //!< queues pour traiter les environnements de façon asynchrone sur le bord du sous-domaine

  // TEST pour amortir le cout des allocs
  BufAddrMng* m_buf_addr_mng=nullptr;

  // Pour synchro algo1
  VarSyncAlgo1* m_vsync_algo1=nullptr;
  Algo1SyncDataMMatDH::PersistentInfo* m_a1_mmat_dh_pi=nullptr;
  Algo1SyncDataMMatD::PersistentInfo* m_a1_mmat_d_pi=nullptr;
};

// Implementation template de computeAndSync
#include "msgpass/ComputeAndSync.h"
#include "msgpass/ComputeMatAndSync.h"
#include "msgpass/EnumerateEnvAndSync.h"

/*---------------------------------------------------------------------------*/
/* Spécialisations pour retourner le nb de fois un type élémentaire DataType */
/* est répété pour un représenter ItemType                                   */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
Integer get_var_degree(MeshVarRefT<ItemType, DataType> var) {
  // Ne devrait jamais être appelé
  ARCANE_ASSERT(false, ("get_var_degree à spécifialiser"));
  return 0;
}
template<typename ItemType, typename DataType>
Integer get_var_degree(MeshVariableScalarRefT<ItemType, DataType> var) {
  return 1;
}

template<typename ItemType, typename DataType>
Integer get_var_degree(MeshVariableArrayRefT<ItemType, DataType> var) {
  return var.arraySize();
}

#endif

