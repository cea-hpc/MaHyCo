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
#include "msgpass/Algo1SyncDataDH.h"
#include "msgpass/Algo1SyncDataD.h"

using namespace Arcane;
using namespace Arcane::Materials;

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

  //! Affecte la version par défaut de implem pour synchronisation var 
  void setDefaultVarSyncVersion(eVarSyncVersion vs_version);

  //! Retourne la version par défaut de implem pour synchronisation var globale
  eVarSyncVersion defaultVarSyncVersion() const;

  //! Buffer d'adresses pour gérer les côuts des allocations
  BufAddrMng* bufAddrMng();

  //! Pour gérer les EnvVarIndex(es) pour les comms
  SyncEnvIndexes* syncEnvIndexes();

  // Retourne l'instance de SyncItems<T> en fonction de T
  template<typename ItemType>
  SyncItems<ItemType>* getSyncItems();

  // Synchronise les éléments fantômes sur une liste de variables
  void synchronize(MeshVariableSynchronizerList& vars, 
    Ref<RunQueue> ref_queue, eVarSyncVersion vs_version=VS_auto);

  // Equivalent à un "var.synchronize()" (implem dépend de vs_version) + plus barrière sur ref_queue
  template<typename MeshVariableRefT>
  void globalSynchronize(Ref<RunQueue> ref_queue, MeshVariableRefT var, eVarSyncVersion vs_version = VS_auto);

  //! Equivalent à un var.synchronize() où var est une variable multi-mat
  template<typename DataType>
  void multiMatSynchronize(CellMaterialVariableScalarRef<DataType> var_menv, Ref<RunQueue> ref_queue, 
      eVarSyncVersion vs_version=VS_auto);

  //! Maj des mailles fantômes d'une liste de variables multi-mat
  void multiMatSynchronize(MeshVariableSynchronizerList& vars, Ref<RunQueue> ref_queue, 
      eVarSyncVersion vs_version=VS_auto);

  /* computeAndSync[OnEvents] */

  // Overlapping entre calcul et communications pour variable multi-mat
  template<typename ItemType, typename Func, typename DataType>
  void computeAndSync(ItemGroupT<ItemType> item_group, Func func, 
      CellMaterialVariableScalarRef<DataType> var, 
      eVarSyncVersion vs_version=VS_auto);

  //! Une fois déroulés les événements de depends_on_evts, 
  // overlapping entre calcul et communications pour variable multi-mat
  template<typename ItemType, typename Func, typename DataType>
  void computeAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    ItemGroupT<ItemType> item_group, Func func, 
    CellMaterialVariableScalarRef<DataType> var, 
    eVarSyncVersion vs_version=VS_auto);

  //! TODO
  template<typename ItemType, typename Func, typename MeshVariableRefT>
  void computeAndSync(ItemGroupT<ItemType> item_group, Func func, 
      MeshVariableRefT var, eVarSyncVersion vs_version=VS_auto);

  // Overlapping entre calcul et communications pour une liste de variables
  template<typename ItemType, typename Func>
  void computeAndSync(ItemGroupT<ItemType> item_group, Func func, 
      MeshVariableSynchronizerList& vars, 
      eVarSyncVersion vs_version=VS_auto);

  //! Une fois déroulés les événements de depends_on_evts, 
  // overlapping entre calcul et communications pour une liste de variables
  template<typename ItemType, typename Func>
  void computeAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
      ItemGroupT<ItemType> item_group, Func func, 
      MeshVariableSynchronizerList& vars, 
      eVarSyncVersion vs_version=VS_auto);

  /* enumerateEnvAndSync[OnEvents] */

  template<typename Func, typename DataType>
  void enumerateEnvAndSync(Func func, CellMaterialVariableScalarRef<DataType> var, 
      eVarSyncVersion vs_version=VS_auto);

  template<typename Func, typename DataType>
  void enumerateEnvAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, CellMaterialVariableScalarRef<DataType> var, 
    eVarSyncVersion vs_version=VS_auto);

  template<typename Func>
  void enumerateEnvAndSync(Func func, MeshVariableSynchronizerList& vars, 
      eVarSyncVersion vs_version=VS_auto);

  template<typename Func>
  void enumerateEnvAndSyncOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    Func func, MeshVariableSynchronizerList& vars, 
    eVarSyncVersion vs_version=VS_auto);

  /* syncAndCompute[OnEvents] */

  // Overlapping entre calcul et communications pour variable multi-mat
  template<typename ItemType, typename Func, typename DataType>
  void syncAndCompute(CellMaterialVariableScalarRef<DataType> var, 
      ItemGroupT<ItemType> item_group, Func func, 
      eVarSyncVersion vs_version=VS_auto);

  //! Une fois déroulés les événements de depends_on_evts, 
  // overlapping entre calcul et communications pour variable multi-mat
  template<typename ItemType, typename Func, typename DataType>
  void syncAndComputeOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
    CellMaterialVariableScalarRef<DataType> var, 
    ItemGroupT<ItemType> item_group, Func func, 
    eVarSyncVersion vs_version=VS_auto);

  //! TODO
  template<typename ItemType, typename Func, typename MeshVariableRefT>
  void syncAndCompute(MeshVariableRefT var, 
      ItemGroupT<ItemType> item_group, Func func, eVarSyncVersion vs_version=VS_auto);

  // Overlapping entre calcul et communications pour une liste de variables
  template<typename ItemType, typename Func>
  void syncAndCompute(MeshVariableSynchronizerList& vars, 
      ItemGroupT<ItemType> item_group, Func func, 
      eVarSyncVersion vs_version=VS_auto);

  //! Une fois déroulés les événements de depends_on_evts, 
  // overlapping entre calcul et communications pour une liste de variables
  template<typename ItemType, typename Func>
  void syncAndComputeOnEvents(ArrayView<Ref<ax::RunQueueEvent>> depends_on_evts,
      MeshVariableSynchronizerList& vars, 
      ItemGroupT<ItemType> item_group, Func func, 
      eVarSyncVersion vs_version=VS_auto);

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
  SyncItems<Face>* m_sync_faces=nullptr;
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

  eVarSyncVersion m_deflt_vs_version;  //!< Version implem pour synchronisation var

  // Pour synchro algo1
  VarSyncAlgo1* m_vsync_algo1=nullptr;
  Algo1SyncDataDH::PersistentInfo* m_a1_dh_pi=nullptr;
  Algo1SyncDataD::PersistentInfo* m_a1_d_pi=nullptr;
};

// Implementation template de computeAndSync
#include "msgpass/ComputeAndSync.h"
#include "msgpass/SyncAndCompute.h"
#include "msgpass/EnumerateEnvAndSync.h"

/*---------------------------------------------------------------------------*/
/* Spécialisations pour retourner le nb de fois un type élémentaire DataType */
/* est répété pour un représenter ItemType                                   */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
Integer get_var_degree(MeshVarRefT<ItemType, DataType> var) {
  // Ne devrait jamais être appelé
  throw NotImplementedException(A_FUNCINFO, String("get_var_degree à spécialiser"));
  return 0;
}
template<typename ItemType, typename DataType>
Integer get_var_degree(MeshVariableScalarRefT<ItemType, DataType> var) {
  (void)var;
  return 1;
}

template<typename ItemType, typename DataType>
Integer get_var_degree(MeshVariableArrayRefT<ItemType, DataType> var) {
  return var.arraySize();
}

#endif

