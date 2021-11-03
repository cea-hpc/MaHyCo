#include "AccEnvDefaultService.h"

#include "cartesian/interface/CellDirectionMng.h"
#include "cartesian/interface/FaceDirectionMng.h"
#include "cartesian/interface/NodeDirectionMng.h"

#include "arcane/materials/CellToAllEnvCellConverter.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/ComponentPartItemVectorView.h"

#include <arcane/AcceleratorRuntimeInitialisationInfo.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

AccEnvDefaultService::AccEnvDefaultService(const ServiceBuildInfo & sbi) : 
  ArcaneAccEnvDefaultObject(sbi), 
  m_node_index_in_cells(platform::getAcceleratorHostMemoryAllocator()),
  m_node_index_in_faces(platform::getAcceleratorHostMemoryAllocator())
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
AccEnvDefaultService::~AccEnvDefaultService() {
  delete m_acc_mem_adv;
  delete m_menv_cell;
  delete m_menv_queue;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
initAcc()
{
  info() << "Using accelerator";
  IApplication* app = subDomain()->application();
  initializeRunner(m_runner,traceMng(),app->acceleratorRuntimeInitialisationInfo());
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
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
_computeNodeIndexInFaces()
{
  info() << "ComputeNodeIndexInFaces";
  // Un noeud est connecté au maximum à max_node_face facettes
  // Calcul pour chaque noeud son index dans chacune des
  // faces à laquelle il est connecté.
  NodeGroup nodes = allNodes();
  Integer nb_node = nodes.size();
  const Integer max_node_face = this->maxNodeFace();
  m_node_index_in_faces.resize(max_node_face*nb_node);
  m_node_index_in_faces.fill(-1);
  auto node_face_cty = this->connectivityView().nodeFace();
  auto face_node_cty = this->connectivityView().faceNode();
  ENUMERATE_NODE(inode,nodes){
    NodeLocalId node = *inode;
    Int32 index = 0;
    Int32 first_pos = node.localId() * max_node_face;
    for( FaceLocalId face : node_face_cty.faces(node) ){
      Int16 node_index_in_face = 0;
      for( NodeLocalId face_node : face_node_cty.nodes(face) ){
        if (face_node==node)
          break;
        ++node_index_in_face;
      }
      m_node_index_in_faces[first_pos + index] = node_index_in_face;
      ++index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/* Pour préparer des données relatives au maillage                           */
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
initMesh(ICartesianMesh* cartesian_mesh)
{
  IMesh* mesh = cartesian_mesh->mesh();

  if (!m_acc_mem_adv) {
    m_acc_mem_adv = new AccMemAdviser(options()->getAccMemAdvise());
  }

  m_connectivity_view.setMesh(mesh);
  // Permet la lecture des cqs quand on boucle sur les noeuds
  _computeNodeIndexInCells();
  _computeNodeIndexInFaces();

  // "Conseils" utilisation de la mémoire unifiée

  m_acc_mem_adv->setReadMostly(m_node_index_in_cells.view());
  m_acc_mem_adv->setReadMostly(m_node_index_in_faces.view());
  
  // CellLocalId
  m_acc_mem_adv->setReadMostly(allCells().view().localIds());
  m_acc_mem_adv->setReadMostly(ownCells().view().localIds());

  // NodeLocalId
  m_acc_mem_adv->setReadMostly(allNodes().view().localIds());
  m_acc_mem_adv->setReadMostly(ownNodes().view().localIds());

  // FaceLocalId
  m_acc_mem_adv->setReadMostly(allFaces().view().localIds());
  m_acc_mem_adv->setReadMostly(ownFaces().view().localIds());

  for(Integer dir(0) ; dir < mesh->dimension() ; ++dir) {
    auto cell_dm = cartesian_mesh->cellDirection(dir);

    // "Conseils" accés mémoire
    m_acc_mem_adv->setReadMostly(cell_dm.innerCells().view().localIds());
    m_acc_mem_adv->setReadMostly(cell_dm.outerCells().view().localIds());
    m_acc_mem_adv->setReadMostly(cell_dm.allCells().view().localIds());
  }
}

/*---------------------------------------------------------------------------*/
/* Calcul des cell_id globaux : permet d'associer à chaque maille impure (mixte) */
/* l'identifiant de la maille globale                                        */
/*---------------------------------------------------------------------------*/
void AccEnvDefaultService::
computeMultiEnvGlobalCellId(IMeshMaterialMng* mesh_material_mng) {
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << "computeMultiEnvGlobalCellId";

  // Calcul des cell_id globaux 
  CellToAllEnvCellConverter all_env_cell_converter(mesh_material_mng);
  ENUMERATE_CELL(icell, allCells()){
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

  m_menv_cell->buildStorage(m_global_cell);

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
}

/*---------------------------------------------------------------------------*/
/* Préparer traitement des environnements sur accélérateur                   */
/*---------------------------------------------------------------------------*/
void AccEnvDefaultService::
initMultiEnv(IMeshMaterialMng* mesh_material_mng) {

  m_menv_queue = new MultiAsyncRunQueue(m_runner, mesh_material_mng->environments().size());

  m_menv_cell = new MultiEnvCellStorage(mesh_material_mng, m_acc_mem_adv);

  updateMultiEnv(mesh_material_mng);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ACCENVDEFAULT(AccEnvDefault, AccEnvDefaultService);

