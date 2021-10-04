#ifndef CARTESIAN_UNSTRUCT_CARTESIAN_MESH_T_H
#define CARTESIAN_UNSTRUCT_CARTESIAN_MESH_T_H

#include "arcane/utils/ArcaneGlobal.h"

#include "cartesian/UnstructConnectivityCellNode.h"
#include "cartesian/UnstructConnectivityCellFaceNode.h"
#include "cartesian/UnstructNeighCells.h"
#include "cartesian/interface/ICartesianMesh.h"
#include "cartesian/interface/CellDirectionMng.h"
#include "cartesian/interface/FaceDirectionMng.h"
#include "cartesian/interface/NodeDirectionMng.h"
#include "cartesian/interface/CartesianConnectivity.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Spécialisation des types pour le non structuré (natif Arcane)
 */
/*---------------------------------------------------------------------------*/
class UnstructCartesianMeshT {
 public:
  // Types qui doivent être définis
  using CellDirectionMngType = CartesianInterface::CellDirectionMng;
  using FaceDirectionMngType = CartesianInterface::FaceDirectionMng;
  using NodeDirectionMngType = CartesianInterface::NodeDirectionMng;
  using CellGroupType = CellGroup;
  using FaceGroupType = FaceGroup;
  using NodeGroupType = NodeGroup;
  using ConnectivityCellNode = Cartesian::UnstructConnectivityCellNode;
  using ConnectivityCellFaceNode = Cartesian::UnstructConnectivityCellFaceNode;
  template<Integer DIM> using NeighCells = Cartesian::UnstructNeighCells<DIM>;
 public:
  UnstructCartesianMeshT(ICartesianMesh* cartesian_mesh) 
  : m_cartesian_mesh (cartesian_mesh),
  m_cc (cartesian_mesh->connectivity()),
  m_trace_mng (cartesian_mesh->mesh()->traceMng()) {
  }
  // Méthodes qui doivent être définies
  //! Gestionnaire de mailles dans la direction idir
  CellDirectionMngType cellDirection(Integer idir) const {
    return m_cartesian_mesh->cellDirection(idir);
  }
  //! Gestionnaire de faces dans la direction idir
  FaceDirectionMngType faceDirection(Integer idir) const {
    return m_cartesian_mesh->faceDirection(idir);
  }
  //! Gestionnaire de noeuds dans la direction idir
  NodeDirectionMngType nodeDirection(Integer idir) const {
    return m_cartesian_mesh->nodeDirection(idir);
  }
  //! Connectivité maille => face => noeuds dans la direction idir
  ConnectivityCellFaceNode connectivityCellFaceNode(Integer idir) {
    return ConnectivityCellFaceNode(m_cartesian_mesh->cellDirection(idir), m_cartesian_mesh->mesh()->dimension());
  }
  //! Connectivité maille => noeuds (sans direction privilégiée)
  ConnectivityCellNode connectivityCellNode() {
    return ConnectivityCellNode(m_cartesian_mesh->mesh());
  }
  //! Connectivité maille => mailles voisines (sans direction privilégiée)
  template<Integer DIM>
  NeighCells<DIM> neighCells() {
    return NeighCells<DIM>(m_cc, m_trace_mng);
  }
 private:
  CartesianInterface::ICartesianMesh* m_cartesian_mesh;
  CartesianInterface::CartesianConnectivity m_cc;
  ITraceMng* m_trace_mng;
};

}

#endif

