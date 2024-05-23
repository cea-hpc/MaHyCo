#include "accenv/MultiEnvMng.h"

#include <arcane/materials/CellToAllEnvCellConverter.h>

/* --------------------------------------------------- */
/* --------------------------------------------------- */
MultiEnvVariableMng::MultiEnvVariableMng(IMeshMaterialMng* mm, ax::Runner& runner, AccMemAdviser* acc_mem_adv) :
  m_mesh_material_mng (mm),
  m_runner (runner),
  m_acc_mem_adv (acc_mem_adv),
  m_max_nb_env (mm->environments().size()),
  
  m_nb_env          (VariableBuildInfo(mm->mesh(), "NbEnv" , IVariable::PNoDump| IVariable::PNoNeedSync)),
  m_l_env_arrays_idx(platform::getAcceleratorHostMemoryAllocator()),
  m_l_env_values_idx(VariableBuildInfo(mm->mesh(), "LEnvValuesIdx" , IVariable::PNoDump| IVariable::PNoNeedSync | IVariable::PSubDomainDepend)),
  m_env_id          (VariableBuildInfo(mm->mesh(), "EnvId" , IVariable::PNoDump| IVariable::PNoNeedSync)),
  m_global_cell     (VariableBuildInfo(mm->mesh(), "GlobalCell" , IVariable::PNoDump| IVariable::PNoNeedSync| IVariable::PSubDomainDepend))
{
  m_l_env_arrays_idx.resize(m_max_nb_env*mm->mesh()->allCells().size());
  acc_mem_adv->setReadMostly(m_l_env_arrays_idx.view());
  m_l_env_values_idx.resize(m_max_nb_env);
}

MultiEnvVariableMng::~MultiEnvVariableMng()
{

}

/*---------------------------------------------------------------------------*/
/* Calcul des cell_id globaux : permet d'associer à chaque maille impure (mixte) */
/* l'identifiant de la maille globale                                        */
/*---------------------------------------------------------------------------*/
void MultiEnvVariableMng::
computeMultiEnvGlobalCellId() 
{
  PROF_ACC_BEGIN(__FUNCTION__);
  m_mesh_material_mng->traceMng()->debug() << "computeMultiEnvGlobalCellId";

  // Calcul des cell_id globaux
  auto async_queue = ax::makeQueueRef(m_runner);
  async_queue->setAsync(true);
  {
    CellToAllEnvCellConverter all_env_cell_converter(m_mesh_material_mng);
    CellToAllEnvCellAccessor cell2allenvcell(m_mesh_material_mng);

    CellGroup all_cells = m_mesh_material_mng->mesh()->allCells();
    auto command = ax::makeCommand(async_queue.get());

    auto out_env_id        = ax::viewOut(command, m_env_id);
    auto out_global_cell   = ax::viewOut(command, m_global_cell);
    auto out_global_cell_g = ax::viewOut(command, m_global_cell.globalVariable());

    command << RUNCOMMAND_ENUMERATE_CELL_ALLENVCELL(cell2allenvcell, cid, all_cells)
    {
      out_global_cell_g[cid] = cid.asInt32();
      Integer nb_env_in_cell = cell2allenvcell.nbEnvironment(cid);
      if (nb_env_in_cell !=1) {
        ENUMERATE_CELL_ALLENVCELL(iev, cid, cell2allenvcell) {
          out_global_cell[*iev] = cid.asInt32();
        }
        // Maille mixte ou vide,
        // Si mixte, contient l'opposé du nombre d'environnements+1
        // Si vide, vaut -1
        out_env_id[cid] = -nb_env_in_cell-1;
      } else {
        AllEnvCell all_env_cell = all_env_cell_converter[cid];
        // Maille pure, cette boucle est de taille 1
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;
          // Cette affectation n'aura lieu qu'une fois
          out_env_id[cid] = ev.environmentId();
        }
      }
    };
  }

  this->asyncBuildStorage(async_queue, m_global_cell);

  asyncCheckMultiEnvGlobalCellId(async_queue);

  async_queue->barrier();

  PROF_ACC_END;
}

void MultiEnvVariableMng::
asyncCheckMultiEnvGlobalCellId([[maybe_unused]] Ref<ax::RunQueue> async_queue) {
#ifdef ARCANE_DEBUG
  m_mesh_material_mng->traceMng()->debug() << "checkMultiEnvGlobalCellId";

  async_queue->barrier();

  // Vérification
  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell(ev.globalCell());
      ARCANE_ASSERT(cell.localId()==m_global_cell[ev], ("lid differents"));
      AllEnvCell all_env_cell(ev.allEnvCell());
      if (all_env_cell.nbEnvironment()==1) {
        ARCANE_ASSERT(m_env_id[cell]==env_id, ("cell pure : environnement id incorrect dans m_env_id[cell]"));
      } else {
        ARCANE_ASSERT(m_env_id[cell]==(-all_env_cell.nbEnvironment()-1), ("cell mixte ou vide : m_env_id[cell] différent de -nbEnvironment()-1"));
      }
    }
  }

  this->asyncCheckStorage(async_queue, m_global_cell);
#endif
}

/*---------------------------------------------------------------------------*/
//! Remplissage de façon asynchrone
/*---------------------------------------------------------------------------*/
void MultiEnvVariableMng::asyncBuildStorage(Ref<ax::RunQueue> async_queue, Materials::MaterialVariableCellInteger& v_global_cell) {
  PROF_ACC_BEGIN(__FUNCTION__);

  Integer max_nb_env = m_max_nb_env; // on ne peut pas utiliser un attribut dans le kernel
  CellToAllEnvCellAccessor cell2allenvcell(m_mesh_material_mng);
  {
    auto command = makeCommand(async_queue.get());

    auto out_nb_env = ax::viewOut(command, m_nb_env);
    auto out_l_env_values_idx = ax::viewOut(command, m_l_env_values_idx);
    auto out_l_env_arrays_idx = m_l_env_arrays_idx.span();

    command << RUNCOMMAND_ENUMERATE_CELL_ALLENVCELL(cell2allenvcell, cid, m_mesh_material_mng->mesh()->allCells()) 
    {
      // Init du nb d'env par maille qui va être calculé à la boucle suivante
      Integer index_cell = 0;

      ENUMERATE_CELL_ALLENVCELL(iev, cid, cell2allenvcell) {
        MatVarIndex mvi = (*iev).localId();

        out_l_env_arrays_idx[cid * max_nb_env + index_cell] = mvi.arrayIndex();
        out_l_env_values_idx[cid][index_cell] = mvi.valueIndex();

        index_cell++;
      }
      out_nb_env[cid] = index_cell;

#ifndef ARCCORE_DEVICE_CODE
      ARCANE_ASSERT(index_cell == cell2allenvcell.nbEnvironment(cid), ("index_cell != cell2allenvcell.nbEnvironment(cid)"));
#endif
    }; // asynchrone, non-bloquant par rapport au CPU
  }

  asyncCheckStorage(async_queue, v_global_cell);
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
//! Verification
/*---------------------------------------------------------------------------*/
void MultiEnvVariableMng::asyncCheckStorage([[maybe_unused]] Ref<ax::RunQueue> async_queue, [[maybe_unused]] Materials::MaterialVariableCellInteger& v_global_cell) {
#ifdef ARCANE_DEBUG
  PROF_ACC_BEGIN(__FUNCTION__);
  async_queue->barrier();
  // m_env_id doit être calculé
  MultiEnvVar<Integer> menv_global_cell(v_global_cell, m_mesh_material_mng);
  auto in_menv_global_cell(menv_global_cell.span());

  ENUMERATE_CELL(icell, m_mesh_material_mng->mesh()->allCells()){
    Cell cell = * icell;
    Integer cid = cell.localId();

    Integer nb_env = m_nb_env[icell];
    Integer nb_env_bis = (m_env_id[cell]>=0 ? 1 : -m_env_id[cell]-1);
    ARCANE_ASSERT(nb_env_bis==nb_env, ("nb_env_bis!=nb_env"));

    for(Integer ienv=0 ; ienv<nb_env ; ++ienv) {
      Integer i=m_l_env_arrays_idx[cid*m_max_nb_env+ienv];
      Integer j=m_l_env_values_idx[icell][ienv];
      EnvVarIndex evi(i,j);

      Integer cid_bis=in_menv_global_cell[evi];
      ARCANE_ASSERT(cid_bis==cid, ("cid_bis!=cid"));
    }
  }
  PROF_ACC_END;
#endif
}


/* --------------------------------------------------- */
/* --------------------------------------------------- */
MultiEnvMng::MultiEnvMng(IMeshMaterialMng* mm, ax::Runner& runner, VarSyncMng* vsync_mng, AccMemAdviser* acc_mem_adv, bool menv_var_mng_needed) :
  m_mesh_material_mng (mm),
  m_runner (runner),
  m_acc_mem_adv (acc_mem_adv),
  m_is_active_cell  (VariableBuildInfo(mm->mesh(), "IsActiveCell" , IVariable::PNoNeedSync)),
  m_is_active_node  (VariableBuildInfo(mm->mesh(), "IsActiveNode" , IVariable::PNoNeedSync))
{
  if (menv_var_mng_needed) {
    m_menv_var_mng = new MultiEnvVariableMng(mm, runner, acc_mem_adv);
  }

  m_menv_queue = new MultiAsyncRunQueue(runner, m_mesh_material_mng->environments().size());

  // 7 = toutes les variables sont synchronisées simultanément 
  m_mesh_material_mng->setSynchronizeVariableVersion(7);
  vsync_mng->initSyncMultiEnv(m_mesh_material_mng);

  updateMultiEnv(vsync_mng);
}

MultiEnvMng::~MultiEnvMng()
{
  delete m_menv_var_mng;
  delete m_menv_queue;
}


void MultiEnvMng::
checkMultiEnvGlobalCellId() {
#ifdef ARCANE_DEBUG
  if (m_menv_var_mng) {
    auto async_queue = ax::makeQueueRef(m_runner);
    m_menv_var_mng->asyncCheckMultiEnvGlobalCellId(async_queue);
    async_queue->barrier();
  }
#endif
}

/*---------------------------------------------------------------------------*/
/* Préparer les données multi-environnement pour l'accélérateur              */
/* A appeler quand la carte des environnements change                        */
/*---------------------------------------------------------------------------*/
void MultiEnvMng::
updateMultiEnv(VarSyncMng* vsync_mng) {
  m_mesh_material_mng->traceMng()->debug() << "updateMultiEnv";

  // Il faut recalculer m_global_cell et m_env_id car la
  // disposition des environnements a changé sur le maillage
  if (m_menv_var_mng)
    m_menv_var_mng->computeMultiEnvGlobalCellId();

  // "Conseils" utilisation de la mémoire unifiée
  ENUMERATE_ENV(ienv,m_mesh_material_mng){
    IMeshEnvironment* env = *ienv;
    m_acc_mem_adv->setReadMostly(env->pureEnvItems().valueIndexes());
    m_acc_mem_adv->setReadMostly(env->impureEnvItems().valueIndexes());
  }

  // Pour mettre à jours des listes pour les comms multi-env
  vsync_mng->updateSyncMultiEnv();
}

/*---------------------------------------------------------------------------*/
// Remplit les tableaux IsActiveCell et IsActiveNode à partir des groupes d'items actifs
/*---------------------------------------------------------------------------*/
/* Chaine de caractères correspondant à un type d'item
 */
template<typename ItemType>
String get_item_string()
{
  return "UNDEF";
}

template<>
String get_item_string<Cell>()
{
  return "cell";
}

template<>
String get_item_string<Node>()
{
  return "node";
}

/* Affecte de façon asynchrone is_active[lid] à true pour les lid \in active_items
 */
template<typename ItemType>
Ref<RunQueue> async_set_active(Runner& runner,
    MeshVariableScalarRefT<ItemType, Byte> is_active,
    ItemGroupT<ItemType> active_items,
    ItemGroupT<ItemType> all_items)
{
  String str_item = get_item_string<ItemType>();

  auto ref_queue = makeQueueRef(runner);
  ref_queue->setAsync(true);
  // ref_queue est détachée par rapport à l'hôte et devra être synchronisée avec barrier()
  // Néanmoins, les 2 kernels doivent s'exécuter l'un après l'autre sur la même queue ref_queue
  // On initialise d'abord à false sur l'ensemble des items (all_items)
  {
    auto command = makeCommand(ref_queue.get());
    auto out_is_active = ax::viewOut(command, is_active);

    // Il s'agit de la macro RUNCOMMAND_ENUMERATE expandée car cette macro ne gère pas bien le type template ItemType
    command.addKernelName(String("init_active_"+str_item)) << A_FUNCINFO 
      << all_items << [=] ARCCORE_HOST_DEVICE (typename ItemType::LocalIdType lid) 
    {
      out_is_active[lid] = false;
    };
  }
  // Puis à true uniquement pour les items dans active_items
  {
    auto command = makeCommand(ref_queue.get());
    auto out_is_active = ax::viewOut(command, is_active);

    command.addKernelName(String("set_active_"+str_item)) << A_FUNCINFO 
      << active_items << [=] ARCCORE_HOST_DEVICE (typename ItemType::LocalIdType lid) 
    {
      out_is_active[lid] = true;
    };
  }

  return ref_queue;
}

void MultiEnvMng::setActiveItemsFromGroups(CellGroup active_cells, NodeGroup active_nodes) 
{
  PROF_ACC_BEGIN(__FUNCTION__);

  m_acc_mem_adv->setReadMostly(active_cells.view().localIds());
  m_acc_mem_adv->setReadMostly(active_nodes.view().localIds());

  // Maj des tableaux m_is_active_cell et m_is_active_node de façon asynchrone
  auto ref_queue_cell = async_set_active<Cell>(m_runner, m_is_active_cell, active_cells, m_mesh_material_mng->mesh()->allCells());
  auto ref_queue_node = async_set_active<Node>(m_runner, m_is_active_node, active_nodes, m_mesh_material_mng->mesh()->allNodes());

  ref_queue_cell->barrier();
  ref_queue_node->barrier();

  PROF_ACC_END;
}

