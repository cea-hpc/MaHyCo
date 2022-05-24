#include "msgpass/VarSyncMng.h"
#include "msgpass/PackTransfer.h"

#include <arcane/IVariableSynchronizer.h>
#include <arcane/IItemFamily.h>
#include <arcane/IParallelMng.h>
#include <arcane/utils/NotSupportedException.h>

// Définie ailleurs
bool is_comm_device_aware();

/*---------------------------------------------------------------------------*/
/* Gère les synchronisations des mailles fantômes par Message Passing        */
/*---------------------------------------------------------------------------*/
VarSyncMng::VarSyncMng(IMesh* mesh, ax::Runner& runner, AccMemAdviser* acc_mem_adv) :
  m_mesh        (mesh),
  m_acc_mem_adv (acc_mem_adv),
  m_runner      (runner)
{
  m_is_device_aware = is_comm_device_aware();

  IItemFamily* cell_family = m_mesh->cellFamily();
  IVariableSynchronizer* var_sync = cell_family->allItemsSynchronizer();

  m_pm = m_mesh->parallelMng();

  // Hypothèse, la liste des voisins est la même quelle que soit le type d'item
  // Donc, je peux récupérer celle issue des mailles
  m_neigh_ranks = var_sync->communicatingRanks();
  m_nb_nei = m_neigh_ranks.size();

  m_sync_cells = new SyncItems<Cell>(m_mesh,m_neigh_ranks, m_acc_mem_adv);
  m_sync_faces = new SyncItems<Face>(m_mesh,m_neigh_ranks, m_acc_mem_adv);
  m_sync_nodes = new SyncItems<Node>(m_mesh,m_neigh_ranks, m_acc_mem_adv);
  m_sync_buffers = new SyncBuffers(isAcceleratorAvailable());
  m_neigh_queues = new MultiAsyncRunQueue(m_runner, m_nb_nei, /*unlimited=*/true);

  _preAllocBuffers();

  m_pack_events.resize(m_nb_nei);
  m_transfer_events.resize(m_nb_nei);
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    m_pack_events[inei]     = makeEventRef(m_runner);
    m_transfer_events[inei] = makeEventRef(m_runner);
  }

  // Queue avec priorité "standard"
  m_ref_queue_inr  = AcceleratorUtils::refQueueAsync(m_runner, QP_default);
  // la priorité doit être la même que celle de la queue qui servira au pack/unpack des buffers de comms = QP_high
  m_ref_queue_bnd  = AcceleratorUtils::refQueueAsync(m_runner, QP_high);
  m_ref_queue_data = AcceleratorUtils::refQueueAsync(m_runner, QP_high);

  // Version par défaut de implem pour synchronisation var globale
  setDefaultGlobVarSyncVersion((isAcceleratorAvailable() ? VS_overlap_evqueue : VS_bulksync_std));

  // Pour synchro algo1
  m_vsync_algo1 = new VarSyncAlgo1(m_pm, m_neigh_ranks);
  m_a1_glob_dh_pi = 
    new Algo1SyncDataGlobDH::PersistentInfo(m_nb_nei, m_runner, m_sync_buffers);
  m_a1_glob_d_pi = 
    new Algo1SyncDataGlobD::PersistentInfo(m_is_device_aware, m_nb_nei, m_runner, m_sync_buffers);
}

VarSyncMng::~VarSyncMng() {
  delete m_sync_cells;
  delete m_sync_faces;
  delete m_sync_nodes;
  delete m_sync_buffers;
  delete m_neigh_queues;

  delete m_sync_evi;

  delete m_menv_queue_inr;
  delete m_menv_queue_bnd;

  delete m_buf_addr_mng;

  delete m_vsync_algo1;
  delete m_a1_mmat_dh_pi;
  delete m_a1_mmat_d_pi;
  delete m_a1_glob_dh_pi;
  delete m_a1_glob_d_pi;
}

/*---------------------------------------------------------------------------*/
/* Initialise les futures synchros multi-mat                                 */
/*---------------------------------------------------------------------------*/
void VarSyncMng::initSyncMultiEnv(IMeshMaterialMng* mesh_material_mng) {
  m_mesh_material_mng = mesh_material_mng;
  if (!m_sync_evi) {
    m_sync_evi = new SyncEnvIndexes(
        MatVarSpace::MaterialAndEnvironment, m_mesh_material_mng,
        m_neigh_ranks, m_sync_cells, m_acc_mem_adv);

    // m_sync_evi doit être créé pour construire m_a1_*
    // _dh_ = Device-Host
    m_a1_mmat_dh_pi = 
      new Algo1SyncDataMMatDH::PersistentInfo(m_nb_nei, m_runner, m_sync_evi, m_sync_buffers);
    // _d_ = only Device
    m_a1_mmat_d_pi = 
      new Algo1SyncDataMMatD::PersistentInfo(m_is_device_aware,
          m_nb_nei, m_runner, m_sync_evi, m_sync_buffers);
  }

  // Pour les traitements multi-env "intérieurs" ou de "bord"
  if (!m_menv_queue_inr) {
    m_menv_queue_inr = new MultiAsyncRunQueue(m_runner, 
        m_mesh_material_mng->environments().size(),
        /*unlimited=*/false, QP_default);
  }
  if (!m_menv_queue_bnd) {
    m_menv_queue_bnd = new MultiAsyncRunQueue(m_runner, 
        m_mesh_material_mng->environments().size(),
        /*unlimited=*/false, QP_high);
  }

  // Buffers mémoire pré-alloués pour minimiser coût des allocs, c'est un TEST
  if (!m_buf_addr_mng) {
    m_buf_addr_mng = new BufAddrMng(m_runner, m_mesh_material_mng);
  }
}

/*---------------------------------------------------------------------------*/
/* Remet à jour les synchros multi-mat quand la carte des environnements a changé */
/*---------------------------------------------------------------------------*/
void VarSyncMng::updateSyncMultiEnv() {
  if (m_sync_evi) {
    m_sync_evi->updateEnvIndexes();
  }
}

/*---------------------------------------------------------------------------*/
/* Retourne vrai si un GPU est dispo pour exécuter les calculs               */
/*---------------------------------------------------------------------------*/
bool VarSyncMng::isAcceleratorAvailable() const {
  return AcceleratorUtils::isAvailable(m_runner);
}

/*---------------------------------------------------------------------------*/
/* Retourne vrai si on peut utiliser les adresses dans DEVICE pour les comms */
/*---------------------------------------------------------------------------*/
bool VarSyncMng::isDeviceAware() const {
  return m_is_device_aware;
}

/*---------------------------------------------------------------------------*/
/* Affecte la version par défaut de implem pour synchronisation var globale  */
/*---------------------------------------------------------------------------*/
void VarSyncMng::setDefaultGlobVarSyncVersion(eVarSyncVersion vs_version) {
  if (vs_version == VS_auto)
  {
    // Détermination automatique
    if (isAcceleratorAvailable())
      m_glob_deflt_vs_version = VS_overlap_evqueue;
    else
      m_glob_deflt_vs_version = VS_bulksync_std; // .synchronize() Arcane
  }
  else
  {
    // Choix imposé par l'utilisateur
    m_glob_deflt_vs_version = vs_version;
  }
}

/*---------------------------------------------------------------------------*/
/* Retourne la version par défaut de implem pour synchronisation var globale */
/*---------------------------------------------------------------------------*/
eVarSyncVersion VarSyncMng::defaultGlobVarSyncVersion() const {
  return m_glob_deflt_vs_version;
}

/*---------------------------------------------------------------------------*/
/* Buffer d'adresses pour gérer les côuts des allocations                    */
/*---------------------------------------------------------------------------*/
BufAddrMng* VarSyncMng::bufAddrMng() {
  return m_buf_addr_mng;
}

/*---------------------------------------------------------------------------*/
/* Pour gérer les EnvVarIndex(es) pour les comms                             */
/*---------------------------------------------------------------------------*/
SyncEnvIndexes* VarSyncMng::syncEnvIndexes() {
  return m_sync_evi;
}

/*---------------------------------------------------------------------------*/
/* Spécialisations pour retourner l'instance de SyncItems<T> en fonction de T*/
/*---------------------------------------------------------------------------*/
template<typename ItemType>
SyncItems<ItemType>* VarSyncMng::getSyncItems() {
  throw NotSupportedException(A_FUNCINFO, "Not implemented for <ItemType>");
  return nullptr;
}

template<>
SyncItems<Cell>* VarSyncMng::getSyncItems() {
  return m_sync_cells;
}

template<>
SyncItems<Face>* VarSyncMng::getSyncItems() {
  return m_sync_faces;
}

template<>
SyncItems<Node>* VarSyncMng::getSyncItems() {
  return m_sync_nodes;
}

/*---------------------------------------------------------------------------*/
/* Effectue une première allocation des buffers pour les communications      */
/* Ceci est une pré-allocation pour miniser le nb de réallocations           */
/*---------------------------------------------------------------------------*/
void VarSyncMng::_preAllocBuffers() {
  auto sync_items = getSyncItems<Cell>();
  
  auto nb_owned_item_idx_pn = sync_items->nbOwnedItemIdxPn();
  auto nb_ghost_item_idx_pn = sync_items->nbGhostItemIdxPn();

  // Alloue pour un buffer de 24 Real3 par maille
  Integer degree = 24;

  m_sync_buffers->resetBuf();
  // On prévoit une taille max du buffer qui va contenir tous les messages
  m_sync_buffers->addEstimatedMaxSz<Real3>(nb_owned_item_idx_pn, degree);
  m_sync_buffers->addEstimatedMaxSz<Real3>(nb_ghost_item_idx_pn, degree);
  // Le buffer de tous les messages est réalloué si pas assez de place
  m_sync_buffers->allocIfNeeded();
}

/*---------------------------------------------------------------------------*/
/* Maj des mailles fantômes d'une liste de variables multi-mat               */
/*---------------------------------------------------------------------------*/
void VarSyncMng::multiMatSynchronize(MeshVariableSynchronizerList& vars, 
    Ref<RunQueue> ref_queue, eVarSyncVersion vs_version)
{
  IAlgo1SyncData* sync_data=nullptr;
  if (vs_version==VS_bulksync_evqueue || vs_version==VS_overlap_evqueue) 
  {
    sync_data = new Algo1SyncDataMMatDH(vars, ref_queue, *m_a1_mmat_dh_pi);
  } 
  else if (vs_version == VS_bulksync_evqueue_d || vs_version==VS_overlap_evqueue_d) 
  {
    sync_data = new Algo1SyncDataMMatD(vars, ref_queue, *m_a1_mmat_d_pi);
  } 
  else 
  {
    throw NotSupportedException(A_FUNCINFO, 
        String::format("Invalid eVarSyncVersion for this method ={0}",(int)vs_version));
  }
  m_vsync_algo1->synchronize(sync_data);
  delete sync_data;
}

/*---------------------------------------------------------------------------*/
/* Maj des mailles fantômes d'une liste de variables (pour l'instant globales) */
/*---------------------------------------------------------------------------*/
void VarSyncMng::synchronize(MeshVariableSynchronizerList& vars, 
    Ref<RunQueue> ref_queue, eVarSyncVersion vs_version)
{
  PROF_ACC_BEGIN(__FUNCTION__);

  if (vs_version == VS_auto) {
    vs_version = defaultGlobVarSyncVersion();
  }

  if (vs_version == VS_bulksync_std)
  {
    // On va construire autant de liste de variables qu'il y a de types d'items
    constexpr Integer MAX_ItemKind = IK_DoF;
    UniqueArray<VariableCollection> all_vc(MAX_ItemKind);

    auto lvars = vars.varsList();
    for(auto var : lvars) {
      if (!var->materialVariable()) {
	// Si on est ici, c'est que la variable est globale
	IVariable* v = var->variable();
	eItemKind item_kind = v->itemKind();
	all_vc[item_kind].add(v);
      }
    }
    // Synchronisation type d'items par type d'item
    for(Integer ik = 0 ; ik<MAX_ItemKind ; ++ik) {
      eItemKind item_kind = static_cast<eItemKind>(ik);
      IItemFamily* item_family = m_mesh->itemFamily(item_kind);
      item_family->synchronize(all_vc[ik]);
    }
  }
  else
  {
    IAlgo1SyncData* sync_data=nullptr;
    if (vs_version==VS_bulksync_evqueue || vs_version==VS_overlap_evqueue) 
    {
      sync_data = new Algo1SyncDataGlobDH(vars, ref_queue, *m_a1_glob_dh_pi);
    } 
    else if (vs_version==VS_bulksync_evqueue_d || vs_version==VS_overlap_evqueue_d) 
    {
      sync_data = new Algo1SyncDataGlobD(vars, ref_queue, *m_a1_glob_d_pi);
    } 
    else 
    {
      throw NotSupportedException(A_FUNCINFO, 
	  String::format("Invalid eVarSyncVersion for this method ={0}",(int)vs_version));
    }
    m_vsync_algo1->synchronize(sync_data);
    delete sync_data;
  }
  
  PROF_ACC_END;
}

