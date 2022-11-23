#include "accenv/MultiEnvMng.h"

#include <arcane/materials/CellToAllEnvCellConverter.h>

MultiEnvMng::MultiEnvMng(IMeshMaterialMng* mm, ax::Runner& runner, VarSyncMng* vsync_mng, AccMemAdviser* acc_mem_adv) :
  m_mesh_material_mng (mm),
  m_runner (runner),
  m_acc_mem_adv (acc_mem_adv),
  m_max_nb_env (mm->environments().size()),

  m_nb_env          (VariableBuildInfo(mm->mesh(), "NbEnv" , IVariable::PNoDump| IVariable::PNoNeedSync)),
  m_l_env_arrays_idx(platform::getAcceleratorHostMemoryAllocator()),
  m_l_env_values_idx(VariableBuildInfo(mm->mesh(), "LEnvValuesIdx" , IVariable::PNoDump| IVariable::PNoNeedSync | IVariable::PSubDomainDepend)),
  m_env_id          (VariableBuildInfo(mm->mesh(), "EnvId" , IVariable::PNoDump| IVariable::PNoNeedSync)),
  m_global_cell     (VariableBuildInfo(mm->mesh(), "GlobalCell" , IVariable::PNoDump| IVariable::PNoNeedSync| IVariable::PSubDomainDepend)),
  m_is_active_cell  (VariableBuildInfo(mm->mesh(), "IsActiveCell" , IVariable::PNoDump| IVariable::PNoNeedSync)),
  m_is_active_node  (VariableBuildInfo(mm->mesh(), "IsActiveNode" , IVariable::PNoDump| IVariable::PNoNeedSync))
{
  m_l_env_arrays_idx.resize(m_max_nb_env*mm->mesh()->allCells().size());
  acc_mem_adv->setReadMostly(m_l_env_arrays_idx.view());
  m_l_env_values_idx.resize(m_max_nb_env);

  m_menv_queue = new MultiAsyncRunQueue(runner, m_mesh_material_mng->environments().size());

  updateMultiEnv(vsync_mng);
  
  // 6 = toutes les variables sont synchronisées simultanément
  m_mesh_material_mng->setSynchronizeVariableVersion(6);
  vsync_mng->initSyncMultiEnv(m_mesh_material_mng);
}

MultiEnvMng::~MultiEnvMng()
{
  delete m_menv_queue;
}

/*---------------------------------------------------------------------------*/
/* Calcul des cell_id globaux : permet d'associer à chaque maille impure (mixte) */
/* l'identifiant de la maille globale                                        */
/*---------------------------------------------------------------------------*/
void MultiEnvMng::
computeMultiEnvGlobalCellId() {
  PROF_ACC_BEGIN(__FUNCTION__);
  m_mesh_material_mng->traceMng()->debug() << "computeMultiEnvGlobalCellId";

  ParallelLoopOptions options;
  options.setPartitioner(ParallelLoopOptions::Partitioner::Auto);

  // Calcul des cell_id globaux 
  arcaneParallelForeach(m_mesh_material_mng->mesh()->allCells(), options, [&](CellVectorView cells) {
    CellToAllEnvCellConverter all_env_cell_converter(m_mesh_material_mng);
    ENUMERATE_CELL(icell, cells)
    {
      Cell cell = * icell;
      Integer cell_id = cell.localId();
      m_global_cell[cell] = cell_id;
      AllEnvCell all_env_cell = all_env_cell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) {
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;
          m_global_cell[ev] = cell_id;
        }
        // Maille mixte ou vide,
        // Si mixte, contient l'opposé du nombre d'environnements+1
        // Si vide, vaut -1
        m_env_id[icell] = -all_env_cell.nbEnvironment()-1;
      } else {
        // Maille pure, cette boucle est de taille 1
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;
          // Cette affectation n'aura lieu qu'une fois
          m_env_id[icell] = ev.environmentId();
        }
      }
    }
  });

  this->buildStorage(m_runner, m_global_cell);

  checkMultiEnvGlobalCellId();
  PROF_ACC_END;
}

void MultiEnvMng::
checkMultiEnvGlobalCellId() {
#ifdef ARCANE_DEBUG
  m_mesh_material_mng->traceMng()->debug() << "checkMultiEnvGlobalCellId";

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

  this->checkStorage(m_global_cell);
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
  computeMultiEnvGlobalCellId();

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
//! Remplissage
/*---------------------------------------------------------------------------*/
void MultiEnvMng::buildStorage(ax::Runner& runner, Materials::MaterialVariableCellInteger& v_global_cell) {
  PROF_ACC_BEGIN(__FUNCTION__);

  auto queue = makeQueue(runner);
  {
    auto command = makeCommand(queue);

    auto inout_nb_env = ax::viewInOut(command, m_nb_env);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, m_mesh_material_mng->mesh()->allCells()){
      // Init du nb d'env par maille qui va être calculé à la boucle suivante
      inout_nb_env[cid] = 0;
    };
  }

  Integer max_nb_env = m_max_nb_env; // on ne peut pas utiliser un attribut dans le kernel
  ENUMERATE_ENV(ienv, m_mesh_material_mng) {
    IMeshEnvironment* env = *ienv;
    Integer env_id = env->id();

    // Mailles mixtes
    {
      auto command = makeCommand(queue);

      auto inout_nb_env = ax::viewInOut(command, m_nb_env);
      auto out_l_env_values_idx = ax::viewOut(command, m_l_env_values_idx);
      auto out_l_env_arrays_idx = m_l_env_arrays_idx.span();

      Span<const Integer> in_global_cell(envView(v_global_cell, env));

      Span<const Int32> in_imp_idx(env->impureEnvItems().valueIndexes());
      Integer nb_imp = in_imp_idx.size();

      command << RUNCOMMAND_LOOP1(iter, nb_imp) {
	auto imix = in_imp_idx[iter()[0]]; // iter()[0] \in [0,nb_imp[
	CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

	Integer index_cell = inout_nb_env[cid];

	// On relève le numéro de l'environnement 
	// et l'indice de la maille dans la liste de mailles mixtes env
	out_l_env_arrays_idx[cid*max_nb_env+index_cell] = env_id+1; // décalage +1 car 0 est pris pour global
	out_l_env_values_idx[cid][index_cell] = imix;

	inout_nb_env[cid] = index_cell+1; // ++ n'est pas supporté
      };
    }

    // Mailles pures
    {
      auto command = makeCommand(queue);

      const auto& pure_env_items = env->pureEnvItems();
      // Pour les mailles pures, valueIndexes() est la liste des ids locaux des mailles
      Span<const Int32> in_cell_id(pure_env_items.valueIndexes());

      auto inout_nb_env = ax::viewInOut(command, m_nb_env);
      auto out_l_env_values_idx = ax::viewOut(command, m_l_env_values_idx);
      auto out_l_env_arrays_idx = m_l_env_arrays_idx.span();

      // Nombre de mailles pures de l'environnement
      Integer nb_pur = pure_env_items.nbItem();

      command << RUNCOMMAND_LOOP1(iter, nb_pur) {
	auto [ipur] = iter(); // ipur \in [0,nb_pur[
	CellLocalId cid(in_cell_id[ipur]); // accés indirect à la valeur de la maille

	Integer index_cell = inout_nb_env[cid];
#ifndef ARCCORE_DEVICE_CODE
	ARCANE_ASSERT(index_cell==0, ("Maille pure mais index_cell!=0"));
#endif
	out_l_env_arrays_idx[cid*max_nb_env+index_cell] = 0; // 0 référence le tableau global
	out_l_env_values_idx[cid][index_cell] = cid.localId();

	// Equivalent à affecter 1
	inout_nb_env[cid] = index_cell+1; // ++ n'est pas supporté
      };
    }
  }

  checkStorage(v_global_cell);
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
//! Verification
/*---------------------------------------------------------------------------*/
void MultiEnvMng::checkStorage([[maybe_unused]] Materials::MaterialVariableCellInteger& v_global_cell) {
#ifdef ARCANE_DEBUG
  PROF_ACC_BEGIN(__FUNCTION__);
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

void MultiEnvMng::setActiveItemsFromGroups(CellGroup active_cells, NodeGroup active_nodes) {

  // Maj des tableaux m_is_active_cell et m_is_active_node de façon asynchrone
  auto ref_queue_cell = async_set_active<Cell>(m_runner, m_is_active_cell, active_cells, m_mesh_material_mng->mesh()->allCells());
  auto ref_queue_node = async_set_active<Node>(m_runner, m_is_active_node, active_nodes, m_mesh_material_mng->mesh()->allNodes());

  ref_queue_cell->barrier();
  ref_queue_node->barrier();
}

