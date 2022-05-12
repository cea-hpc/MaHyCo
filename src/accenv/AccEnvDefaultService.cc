#include "AccEnvDefaultService.h"

#include "cartesian/interface/CellDirectionMng.h"
#include "cartesian/interface/FaceDirectionMng.h"
#include "cartesian/interface/NodeDirectionMng.h"

#include "arcane/materials/CellToAllEnvCellConverter.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/ComponentPartItemVectorView.h"

#include <arcane/AcceleratorRuntimeInitialisationInfo.h>
#include <arcane/IParallelMng.h>
#include <arcane/ParallelMngUtils.h>
#include <arcane/IParallelTopology.h>

#if defined(ACCENV_HWLOC) && defined(ARCANE_COMPILING_CUDA)
#warning "HWLOC présent pour placement GPUs"
// Code issu de https://www.open-mpi.org/faq/?category=runcuda
#include <cuda.h>
#include <mpi.h>
#include <hwloc.h>

#define ABORT_ON_ERROR(func)                          \
  { CUresult res;                                     \
    res = func;                                       \
    if (CUDA_SUCCESS != res) {                        \
        printf("%s returned error=%d\n", #func, res); \
        abort();                                      \
    }                                                 \
  }

#define ABORT_ON_HWLOC_ERROR(func)                    \
  { int res;                                          \
    res = func;                                       \
    if (0 != res) {                                   \
        printf("%s returned error=%d\n", #func, res); \
        abort();                                      \
    }                                                 \
  }

/**
 * This function searches for all the GPUs that are hanging off a NUMA
 * node.  It walks through each of the PCI devices and looks for ones
 * with the NVIDIA vendor ID.  It then stores them into an array.
 * Note that there can be more than one GPU on the NUMA node.
 */
 
void find_gpus(hwloc_topology_t topology, hwloc_obj_t parent, hwloc_obj_t child,
    int* gpuIndex, hwloc_obj_t gpus[]) {
  hwloc_obj_t pcidev;
  pcidev = hwloc_get_next_child(topology, parent, child);
  if (NULL == pcidev) {
    return;
  } else if (0 != pcidev->arity) {
    /* This device has children so need to look recursively at them */
    find_gpus(topology, pcidev, NULL, gpuIndex, gpus);
    find_gpus(topology, parent, pcidev, gpuIndex, gpus);
  } else {
    if (pcidev->attr->pcidev.vendor_id == 0x10de) {
      gpus[(*gpuIndex)++] = pcidev;
    }
    find_gpus(topology, parent, pcidev, gpuIndex, gpus);
  }
}
#endif

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

AccEnvDefaultService::AccEnvDefaultService(const ServiceBuildInfo & sbi) : 
  ArcaneAccEnvDefaultObject(sbi), 
  m_node_index_in_cells(platform::getAcceleratorHostMemoryAllocator())
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
AccEnvDefaultService::~AccEnvDefaultService() {
  delete m_acc_mem_adv;
  delete m_menv_cell;
  delete m_menv_queue;
  delete m_vsync_mng;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
initAcc()
{
  info() << "Using accelerator";
  IApplication* app = subDomain()->application();
  if (options()->getHeterogPartition() == HP_none) 
  {
    initializeRunner(m_runner,traceMng(),app->acceleratorRuntimeInitialisationInfo());
  }
  else if (options()->getHeterogPartition() == HP_heterog1)
  {
    Integer rank = mesh()->parallelMng()->commRank();
    bool is_host_only = (rank > 0);
    if (is_host_only) {
      std::ostringstream ostr;
      ostr << "(Host) P=" << mesh()->parallelMng()->commRank();
      if (TaskFactory::isActive()) {
        ostr << " using Task runtime";
        m_runner.initialize(ax::eExecutionPolicy::Thread);
      } else {
        ostr << " using Sequential runtime";
        m_runner.initialize(ax::eExecutionPolicy::Sequential);
      }
      pinfo() << ostr.str();
    } else {
      // Initialisation classique : accelerateur ou pas
      initializeRunner(m_runner,traceMng(),app->acceleratorRuntimeInitialisationInfo());
    }
  }

  bool is_acc_av = AcceleratorUtils::isAvailable(m_runner);
  if (is_acc_av && options()->getDeviceAffinity() == DA_world_rank) 
  {
    info() << "Placement GPU : device =  world_rank%device_count";

    Integer device_count = AcceleratorUtils::deviceCount();
    if (device_count>0) {
      Integer rank = mesh()->parallelMng()->commRank();
      Integer device=rank%device_count;
      AcceleratorUtils::setDevice(device);
      pinfo() << "Processus " << rank  
        << " : Device " << device << " (pour " << device_count << " device(s))";
    }
  }
  if (options()->getDeviceAffinity() == DA_node_rank) 
  {
    IParallelMng* pm = mesh()->parallelMng();
    // Attention, createTopology() est une opération collective
    Ref<IParallelTopology> pt = ParallelMngUtils::createTopologyRef(pm);
    auto node_ranks = pt->machineRanks();
    // Attention, createSubParallelMng() est une opération collective
    Ref<IParallelMng> pm_node = pm->createSubParallelMngRef(node_ranks);

    if (is_acc_av) {
      info() << "Placement GPU : device =  node_rank%device_count";

      Integer device_count = AcceleratorUtils::deviceCount();
      if (device_count>0) {
        Integer rank = pm->commRank();
        Integer node_rank = pm_node->commRank();
        Integer device=node_rank%device_count;
        AcceleratorUtils::setDevice(device);
        pinfo() << "Processus " << rank  << " (node_rank=" << node_rank << ")"
          << " : Device " << device << " (pour " << device_count << " device(s))";
      }
    }
  }
#ifdef ARCANE_COMPILING_CUDA
#ifdef ACCENV_HWLOC
  if (is_acc_av && options()->getDeviceAffinity() == DA_cu_hwloc) 
  {
    info() << "Placement GPU avec hwloc";
    // Code issu de https://www.open-mpi.org/faq/?category=runcuda
    //
#if 0
    const unsigned long flags = HWLOC_TOPOLOGY_FLAG_IO_DEVICES | HWLOC_TOPOLOGY_FLAG_IO_BRIDGES;
#endif
    hwloc_cpuset_t newset;
    hwloc_obj_t node, bridge;
    char pciBusId[16];
    CUdevice dev;
    char devName[256];
    hwloc_topology_t topology = NULL;
    int gpuIndex = 0;
    hwloc_obj_t gpus[16] = {0};

    /* Now decide which GPU to pick.  This requires hwloc to work properly.
     * We first see which CPU we are bound to, then try and find a GPU nearby.
     */
    ABORT_ON_HWLOC_ERROR( hwloc_topology_init(&topology) );
#if 0
    ABORT_ON_HWLOC_ERROR( hwloc_topology_set_flags(topology, flags) );
#else
    ABORT_ON_HWLOC_ERROR( hwloc_topology_set_io_types_filter(topology, HWLOC_TYPE_FILTER_KEEP_IMPORTANT) );
#endif
    ABORT_ON_HWLOC_ERROR( hwloc_topology_load(topology) );
    newset = hwloc_bitmap_alloc();
    ABORT_ON_HWLOC_ERROR( hwloc_get_last_cpu_location(topology, newset, 0) );

    /* Get the object that contains the cpuset */
    node = hwloc_get_first_largest_obj_inside_cpuset(topology, newset);

    /* Climb up from that object until we find the HWLOC_OBJ_NODE */
    while (node->type != HWLOC_OBJ_NODE) {
      node = node->parent;
    }

    /* Now look for the HWLOC_OBJ_BRIDGE.  All PCI busses hanging off the
     * node will have one of these */
    bridge = hwloc_get_next_child(topology, node, NULL);
    while (bridge->type != HWLOC_OBJ_BRIDGE) {
      bridge = hwloc_get_next_child(topology, node, bridge);
    }

    /* Now find all the GPUs on this NUMA node and put them into an array */
    find_gpus(topology, bridge, NULL, &gpuIndex, gpus);

    /* Now select the first GPU that we find */
    if (gpus[0] == 0) {
      printf("No GPU found\n");
      exit(1);
    } else {
      sprintf(pciBusId, "%.2x:%.2x:%.2x.%x", gpus[0]->attr->pcidev.domain, gpus[0]->attr->pcidev.bus,
          gpus[0]->attr->pcidev.dev, gpus[0]->attr->pcidev.func);
      ABORT_ON_ERROR(cuDeviceGetByPCIBusId(&dev, pciBusId));
      ABORT_ON_ERROR(cuDeviceGetName(devName, 256, dev));

      // Affichage
      Integer device, device_count;
      cudaGetDevice(&device_count);
      cudaGetDevice(&device);
      pinfo() << "Processus " << mesh()->parallelMng()->commRank()
        << " : Device " << device << " (pour " << device_count << " device(s))"
        << " (Selected GPU=" << pciBusId << ", name=" << devName << ")";
    }
  }
#endif // ACCENV_HWLOC
#endif // ARCANE_COMPILING_CUDA
}

/*---------------------------------------------------------------------------*/
/* Référence sur une queue asynchrone créée avec un niveau de priorité       */
/*---------------------------------------------------------------------------*/
Ref<ax::RunQueue> AccEnvDefaultService::
refQueueAsync(eQueuePriority qp) {
  return AcceleratorUtils::refQueueAsync(m_runner, qp);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void AccEnvDefaultService::
_computeNodeIndexInCells() {
  debug() << "_computeNodeIndexInCells";
  // Un noeud est connecté au maximum à max_node_cell mailles
  // Calcul pour chaque noeud son index dans chacune des
  // mailles à laquelle il est connecté.
  NodeGroup nodes = allNodes();
  Integer nb_node = nodes.size();
  const Integer max_node_cell = this->maxNodeCell();
  m_node_index_in_cells.resize(max_node_cell*nb_node);
  m_node_index_in_cells.fill(-1);
  auto node_cell_cty = this->connectivityView().nodeCell();
  auto cell_node_cty = this->connectivityView().cellNode();
  ENUMERATE_NODE(inode,nodes){
    NodeLocalId node = *inode;
    Int32 index = 0; 
    Int32 first_pos = node.localId() * max_node_cell;
    for( CellLocalId cell : node_cell_cty.cells(node) ){
      Int16 node_index_in_cell = 0; 
      for( NodeLocalId cell_node : cell_node_cty.nodes(cell) ){
        if (cell_node==node)
          break;
        ++node_index_in_cell;
      }    
      m_node_index_in_cells[first_pos + index] = node_index_in_cell;
      ++index;
    }    
  }
}

/*---------------------------------------------------------------------------*/
/* Pour préparer des données relatives au maillage                           */
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
initMesh(IMesh* mesh)
{
  if (!m_acc_mem_adv) {
    m_acc_mem_adv = new AccMemAdviser(options()->getAccMemAdvise());
  }

  m_connectivity_view.setMesh(mesh);
  // Permet la lecture des cqs quand on boucle sur les noeuds
  _computeNodeIndexInCells();

  // "Conseils" utilisation de la mémoire unifiée

  m_acc_mem_adv->setReadMostly(m_node_index_in_cells.view());
  
  // CellLocalId
  m_acc_mem_adv->setReadMostly(allCells().view().localIds());
  m_acc_mem_adv->setReadMostly(ownCells().view().localIds());

  // NodeLocalId
  m_acc_mem_adv->setReadMostly(allNodes().view().localIds());
  m_acc_mem_adv->setReadMostly(ownNodes().view().localIds());

  // FaceLocalId
  m_acc_mem_adv->setReadMostly(allFaces().view().localIds());
  m_acc_mem_adv->setReadMostly(ownFaces().view().localIds());

  m_vsync_mng = new VarSyncMng(mesh, m_runner, m_acc_mem_adv);
}

/*---------------------------------------------------------------------------*/
/* Calcul des cell_id globaux : permet d'associer à chaque maille impure (mixte) */
/* l'identifiant de la maille globale                                        */
/*---------------------------------------------------------------------------*/
void AccEnvDefaultService::
computeMultiEnvGlobalCellId(IMeshMaterialMng* mesh_material_mng) {
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "computeMultiEnvGlobalCellId";

  ParallelLoopOptions options;
  options.setPartitioner(ParallelLoopOptions::Partitioner::Auto);

  // Calcul des cell_id globaux 
  arcaneParallelForeach(allCells(), options, [&](CellVectorView cells) {
    CellToAllEnvCellConverter all_env_cell_converter(mesh_material_mng);
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

  m_menv_cell->buildStorage(m_runner, m_global_cell);

  checkMultiEnvGlobalCellId(mesh_material_mng);
  PROF_ACC_END;
}

void AccEnvDefaultService::
checkMultiEnvGlobalCellId(IMeshMaterialMng* mesh_material_mng) {
#ifdef ARCANE_DEBUG
  debug() << "checkMultiEnvGlobalCellId";

  // Vérification
  ENUMERATE_ENV(ienv, mesh_material_mng) {
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

  m_menv_cell->checkStorage(m_global_cell);
#endif
}

/*---------------------------------------------------------------------------*/
/* Préparer les données multi-envronnement pour l'accélérateur               */
/* A appeler quand la carte des environnements change                        */
/*---------------------------------------------------------------------------*/
void AccEnvDefaultService::
updateMultiEnv(IMeshMaterialMng* mesh_material_mng) {
  debug() << "updateMultiEnv";

  // Il faut recalculer m_global_cell et m_env_id car la
  // disposition des environnements a changé sur le maillage
  computeMultiEnvGlobalCellId(mesh_material_mng);

  // "Conseils" utilisation de la mémoire unifiée
  ENUMERATE_ENV(ienv,mesh_material_mng){
    IMeshEnvironment* env = *ienv;
    m_acc_mem_adv->setReadMostly(env->pureEnvItems().valueIndexes());
  }

  // Pour mettre à jours des listes pour les comms multi-env
  m_vsync_mng->updateSyncMultiEnv();
}

/*---------------------------------------------------------------------------*/
/* Préparer traitement des environnements sur accélérateur                   */
/*---------------------------------------------------------------------------*/
void AccEnvDefaultService::
initMultiEnv(IMeshMaterialMng* mesh_material_mng) {

  // 6 = toutes les variables sont synchronisées simultanément
  mesh_material_mng->setSynchronizeVariableVersion(6);
  m_vsync_mng->initSyncMultiEnv(mesh_material_mng);

  m_menv_queue = new MultiAsyncRunQueue(m_runner, mesh_material_mng->environments().size());

  m_menv_cell = new MultiEnvCellStorage(mesh_material_mng, m_acc_mem_adv);

  updateMultiEnv(mesh_material_mng);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ACCENVDEFAULT(AccEnvDefault, AccEnvDefaultService);

