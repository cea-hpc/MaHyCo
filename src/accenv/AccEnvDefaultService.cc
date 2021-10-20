#include "AccEnvDefaultService.h"

#include "cartesian/interface/CellDirectionMng.h"
#include "cartesian/interface/FaceDirectionMng.h"
#include "cartesian/interface/NodeDirectionMng.h"

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

  m_connectivity_view.setMesh(mesh);
  // Permet la lecture des cqs quand on boucle sur les noeuds
  _computeNodeIndexInCells();
  _computeNodeIndexInFaces();

#ifdef ARCANE_COMPILING_CUDA
  // "Conseils" utilisation de la mémoire unifiée
  int device = -1;
  cudaGetDevice(&device);

  mem_adv_set_read_mostly(m_node_index_in_cells.view(), device);
  mem_adv_set_read_mostly(m_node_index_in_faces.view(), device);
  
  // CellLocalId
  mem_adv_set_read_mostly(allCells().view().localIds(), device);
  mem_adv_set_read_mostly(ownCells().view().localIds(), device);

  // NodeLocalId
  mem_adv_set_read_mostly(allNodes().view().localIds(), device);
  mem_adv_set_read_mostly(ownNodes().view().localIds(), device);

  // FaceLocalId
  mem_adv_set_read_mostly(allFaces().view().localIds(), device);
  mem_adv_set_read_mostly(ownFaces().view().localIds(), device);

  for(Integer dir(0) ; dir < mesh->dimension() ; ++dir) {
    auto cell_dm = cartesian_mesh->cellDirection(dir);

    // "Conseils" accés mémoire
    mem_adv_set_read_mostly(cell_dm.innerCells().view().localIds(), device);
    mem_adv_set_read_mostly(cell_dm.outerCells().view().localIds(), device);
    mem_adv_set_read_mostly(cell_dm.allCells().view().localIds(), device);
  }
#endif
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ACCENVDEFAULT(AccEnvDefault, AccEnvDefaultService);

