#include "AccEnvDefaultService.h"

#include <arcane/AcceleratorRuntimeInitialisationInfo.h>
#include <arcane/IParallelMng.h>
#include <arcane/ParallelMngUtils.h>
#include <arcane/IParallelTopology.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>

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
  delete m_menv_mng;
  delete m_vsync_mng;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
initAcc()
{
  info() << "Using accelerator";
  if (options()->getHeterogPartition() == HP_none) 
  {
    m_runner_ptr = subDomain()->acceleratorMng()->defaultRunner();
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
        m_host_runner.initialize(ax::eExecutionPolicy::Thread);
      } else {
        ostr << " using Sequential runtime";
        m_host_runner.initialize(ax::eExecutionPolicy::Sequential);
      }
      pinfo() << ostr.str();
      m_runner_ptr = &m_host_runner;
    } else {
      // Initialisation classique : accelerateur ou pas
      m_runner_ptr = subDomain()->acceleratorMng()->defaultRunner();
    }
  }
  else
  {
    m_runner_ptr = subDomain()->acceleratorMng()->defaultRunner();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
initPAcc()
{
  m_pacc_needed = true;
}

/*---------------------------------------------------------------------------*/
/* Référence sur une queue asynchrone créée avec un niveau de priorité       */
/*---------------------------------------------------------------------------*/
Ref<ax::RunQueue> AccEnvDefaultService::
refQueueAsync(eQueuePriority qp) {
  return AcceleratorUtils::refQueueAsync(runner(), qp);
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

  m_vsync_mng = new VarSyncMng(mesh, runner(), m_acc_mem_adv, m_pacc_needed);
  m_vsync_mng->setDefaultVarSyncVersion(options()->getVarSyncVersion());
}

/*---------------------------------------------------------------------------*/
/* Fabrique l'instance de MultiEnvMng                                        */
/*---------------------------------------------------------------------------*/
void AccEnvDefaultService::
createMultiEnvMng(IMeshMaterialMng* mesh_material_mng)
{
  if (m_menv_mng)
  {
    throw FatalErrorException(A_FUNCINFO, "Une instance de MultiEnvMng existe déjà");
  }
  m_menv_mng = new MultiEnvMng(mesh_material_mng, this->runner(), m_vsync_mng, m_acc_mem_adv, m_pacc_needed);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ACCENVDEFAULT(AccEnvDefault, AccEnvDefaultService);

