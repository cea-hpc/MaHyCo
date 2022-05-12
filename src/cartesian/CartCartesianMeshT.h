#ifndef CARTESIAN_CART_CARTESIAN_MESH_T_H
#define CARTESIAN_CART_CARTESIAN_MESH_T_H

#include "arcane/utils/ArcaneGlobal.h"

#include "cartesian/FactCartDirectionMng.h"
#include "cartesian/CartConnectivityCellNode.h"
#include "cartesian/CartConnectivityCellFaceNode.h"
#include "cartesian/CartNeighCells.h"
#include "cartesian/interface/ICartesianMesh.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Spécialisation des types pour le Cartésien
 */
/*---------------------------------------------------------------------------*/
class CartCartesianMeshT {
 public:
  // Types qui doivent être définis
  using CellDirectionMngType = Cartesian::CartCellDirectionMng;
  using FaceDirectionMngType = Cartesian::CartFaceDirectionMng;
  using NodeDirectionMngType = Cartesian::CartNodeDirectionMng;
  using CellGroupType = Cartesian::CartCellGroup;
  using FaceGroupType = Cartesian::CartFaceGroup;
  using NodeGroupType = Cartesian::CartNodeGroup;
  using ConnectivityCellNode = Cartesian::CartConnectivityCellNode;
  using ConnectivityCellFaceNode = Cartesian::CartConnectivityCellFaceNode;
  template<Integer DIM>  using NeighCells = Cartesian::CartNeighCells<DIM>;
 private:
  using CartesianGrid = typename ConnectivityCellFaceNode::CartesianGrid;
 public:
  CartCartesianMeshT(ICartesianMesh* cartesian_mesh) 
  :
  m_fact_cart_dm (cartesian_mesh->mesh()),
  m_cart_grid (*m_fact_cart_dm.cartesianGrid()),
  m_trace_mng (cartesian_mesh->mesh()->traceMng()) {
  }
  // Méthodes qui doivent être définies
  //! Gestionnaire de mailles dans la direction idir
  CellDirectionMngType cellDirection(Integer idir) const {
    return m_fact_cart_dm.cellDirection(idir);
  }
  //! Gestionnaire de faces dans la direction idir
  FaceDirectionMngType faceDirection(Integer idir) const {
    return m_fact_cart_dm.faceDirection(idir);
  }
  //! Gestionnaire de noeuds dans la direction idir
  NodeDirectionMngType nodeDirection(Integer idir) const {
    return m_fact_cart_dm.nodeDirection(idir);
  }
  //! Connectivité maille => face => noeuds dans la direction idir
  ConnectivityCellFaceNode connectivityCellFaceNode(Integer idir) {
    return ConnectivityCellFaceNode(idir, m_cart_grid);
  }
  //! Connectivité maille => noeuds (sans direction privilégiée)
  ConnectivityCellNode connectivityCellNode() {
    return ConnectivityCellNode(m_cart_grid, /*sort_node=*/Cartesian::CartConnectivityCellNode::SN_trigo);
  }
  //! Connectivité maille => mailles voisines (sans direction privilégiée)
  template<Integer DIM>
  NeighCells<DIM> neighCells() {
    return NeighCells<DIM>(m_cart_grid.cartNumCell(), m_trace_mng);
  }
 private:
  Cartesian::FactCartDirectionMng m_fact_cart_dm;
  const CartesianGrid& m_cart_grid;
  ITraceMng* m_trace_mng;
};

}

#endif

