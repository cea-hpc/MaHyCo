#ifndef CARTESIAN_UNSTRUCT_CONNECTIVITY_CELL_FACE_NODE_H
#define CARTESIAN_UNSTRUCT_CONNECTIVITY_CELL_FACE_NODE_H

#include "arcane/Item.h"
#include "arcane/ItemEnumerator.h"
#include "cartesian/interface/CellDirectionMng.h"
#include "cartesian/CartTypes.h"

namespace Cartesian {
  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Implementation connectivité maille => noeud sur face dans dir 
 * à partir d'une connectivité non structurée
 */
/*---------------------------------------------------------------------------*/
class UnstructConnectivityCellFaceNode {
 public:
  using NodeType = NodeLocalId;
  using CellEnumeratorType = ItemEnumeratorT<Cell>;
  
  UnstructConnectivityCellFaceNode(CartesianInterface::CellDirectionMng cell_dm, Integer dim) 
  : m_cell_dm (cell_dm) {
    m_nb_nodes_face_dir = 1 << (dim-1); //! Nb de noeuds sur la face orthogonale a dir

  }

  // Nb de noeuds sur la face orthogonale a dir
  Integer nbNode() const {
    return m_nb_nodes_face_dir;
  }

  void initCartCell(const CellEnumeratorType &c) {
    DirCellFace cf(m_cell_dm.cellFace(*c));
    m_side_face = (m_side == MS_previous ? cf.previous() : cf.next());
  }

  void initSide(eMeshSide side) {
    m_side = side;
  }

  NodeType node(Integer inode) const {
    return m_side_face.node(inode);
  }

 private:
  

 private:
  Integer m_nb_nodes_face_dir; //! Nb de noeuds sur la face orthogonale a dir

  eMeshSide m_side = MS_invalid; //! MS_previous ou MS_next

  Face m_side_face;
  CartesianInterface::CellDirectionMng m_cell_dm;
};

}

#endif

