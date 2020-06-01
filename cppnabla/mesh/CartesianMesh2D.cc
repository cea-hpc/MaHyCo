/*******************************************************************************
 * Copyright (c) 2020 CEA
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
	 * Contributors: see AUTHORS file
 *******************************************************************************/
#include "mesh/CartesianMesh2D.h"

#define __CARTESIAN__
//#define __VALIDATION__


namespace nablalib
{

#ifdef __VALIDATION__
bool areVectorEqual(const vector<int>& a, vector<int>& b)
{
  vector<int> aa(a);
  vector<int> bb(b);
  std::sort(aa.begin(), aa.end());
  std::sort(bb.begin(), bb.end());
  return aa==bb;
}
#endif
  
inline int
CartesianMesh2D::ij2kCell(const int& i, const int&j) const noexcept
{
  return i*m_nb_x_quads + j;
}
  
inline void
CartesianMesh2D::k2ijCell(const int& k, int& i, int&j) const noexcept
{
  i = k / m_nb_x_quads;
  j = k - i*m_nb_x_quads;
}

inline int
CartesianMesh2D::ij2kNode(const int& i, const int&j) const noexcept
{
  return i*(m_nb_x_quads+1) + j;
}
  
inline void
CartesianMesh2D::k2ijNode(const int& k, int& i, int&j) const noexcept
{
  i = k / (m_nb_x_quads+1);
  j = k - i*(m_nb_x_quads+1);
}


  
CartesianMesh2D::CartesianMesh2D(
	MeshGeometry<2>* geometry,
	const vector<int>& inner_nodes_ids,
	const vector<int>& top_nodes_ids,
	const vector<int>& bottom_nodes_ids,
	const vector<int>& left_nodes_ids,
	const vector<int>& right_nodes_ids,
	const int top_left_node_id,
	const int top_right_node_id,
	const int bottom_left_node_id,
	const int bottom_right_node_id)
: m_geometry(geometry)
, m_inner_nodes(inner_nodes_ids)
, m_top_nodes(top_nodes_ids)
, m_bottom_nodes(bottom_nodes_ids)
, m_left_nodes(left_nodes_ids)
, m_right_nodes(right_nodes_ids)
, m_top_left_node(top_left_node_id)
, m_top_right_node(top_right_node_id)
, m_bottom_left_node(bottom_left_node_id)
, m_bottom_right_node(bottom_right_node_id)
, m_nb_x_quads(bottom_nodes_ids.size()+1)
, m_nb_y_quads(left_nodes_ids.size()+1)
{
	// outer faces
	auto edges = m_geometry->getEdges();
	for (int edgeId(0); edgeId < edges.size(); ++edgeId) {
	        m_faces.emplace_back(edgeId);
		if (!isInnerEdge(edges[edgeId]))
			m_outer_faces.emplace_back(edgeId);
		else {
		        m_inner_faces.emplace_back(edgeId);
			if (isInnerVerticalEdge(edges[edgeId])) {
				m_inner_vertical_faces.emplace_back(edgeId);
			}
			else if (isInnerHorizontalEdge(edges[edgeId])) {
				m_inner_horizontal_faces.emplace_back(edgeId);
			}
			else {
				cout << "The inner edge should be either vertical or horizontal \n";
				terminate();
			}
		}
	}
}

const array<int, 4>&
CartesianMesh2D::getNodesOfCell(const int& cellId) const noexcept
{
	return m_geometry->getQuads()[cellId].getNodeIds();
}

const array<int, 2>&
CartesianMesh2D::getNodesOfFace(const int& faceId) const noexcept
{
	return m_geometry->getEdges()[faceId].getNodeIds();
}

int
CartesianMesh2D::getFirstNodeOfFace(const int& faceId) const noexcept
{
	return m_geometry->getEdges()[faceId].getNodeIds()[0];
}
int
CartesianMesh2D::getSecondNodeOfFace(const int& faceId) const noexcept
{
	return m_geometry->getEdges()[faceId].getNodeIds()[1];
}
vector<int>
CartesianMesh2D::getCellsOfNode(const int& nodeId) const noexcept
{
#if defined(__VALIDATION__) || !defined(__CARTESIAN__)
	vector<int> candidateQuadIds;
	auto quads = m_geometry->getQuads();
	for (int quadId(0); quadId < quads.size(); ++quadId) {
		if (find(quads[quadId].getNodeIds().begin(), quads[quadId].getNodeIds().end(), nodeId) != quads[quadId].getNodeIds().end())
			candidateQuadIds.emplace_back(quadId);
	}
#endif

#if defined(__VALIDATION__) || defined(__CARTESIAN__)	
	std::vector<int> cells;
	int i,j;
	k2ijNode(nodeId, i,j);
	if (i < m_nb_y_quads && j < m_nb_x_quads) cells.emplace_back(ij2kCell(i,j));
	if (i < m_nb_y_quads && j > 0)              cells.emplace_back(ij2kCell(i,j-1));
	if (i > 0              && j < m_nb_x_quads) cells.emplace_back(ij2kCell(i-1,j));
	if (i > 0              && j > 0)              cells.emplace_back(ij2kCell(i-1,j-1));
#endif

#ifdef __VALIDATION__
	if (!areVectorEqual(candidateQuadIds, cells))
	{
	  cout << "getCellsOfNode : Nabla and index versions do not match !!\n";
	  std::terminate();
	}
#endif

#ifdef __CARTESIAN__
	return cells;
#else
	return candidateQuadIds;
#endif	
}

vector<int>
CartesianMesh2D::getCellsOfFace(const int& faceId) const
{
#if defined(__VALIDATION__) || !defined(__CARTESIAN__)
	std::vector<int> cellsOfFace;
	const auto& nodes(getNodesOfFace(faceId));
	for (auto nodeId : nodes)
	{
		auto adjacentCells(getCellsOfNode(nodeId));
		for (int quadId : adjacentCells)
			if (getNbCommonIds(nodes, m_geometry->getQuads()[quadId].getNodeIds()) == 2)
				cellsOfFace.emplace_back(quadId);
	}
	std::sort(cellsOfFace.begin(), cellsOfFace.end());
	cellsOfFace.erase(std::unique(cellsOfFace.begin(), cellsOfFace.end()), cellsOfFace.end());
#endif

#if defined(__VALIDATION__) || defined(__CARTESIAN__)
	std::vector<int> cells;	
	int i_f = faceId / (2*m_nb_x_quads+1);
	int k_f = faceId - i_f*(2*m_nb_x_quads+1);
	
	if (i_f < m_nb_y_quads)  // all except upper bound faces
	{
	  if (k_f == 2*m_nb_x_quads) // right bound edge 
	  {
	    cells.emplace_back(ij2kCell(i_f, m_nb_x_quads-1));
	  }
          else if (k_f == 1) // left bound edge
	  {
	    cells.emplace_back(ij2kCell(i_f, 0));
	  }
	  else if (k_f % 2 == 0) // horizontal edge
	  {
	    cells.emplace_back(ij2kCell(i_f, k_f/2));
	    if (i_f > 0) // Not bottom bound edge
	      cells.emplace_back(ij2kCell(i_f-1, k_f/2));
	  } else // vertical edge (neither left bound nor right bound)
	  {
	    cells.emplace_back(ij2kCell(i_f, (k_f-1)/2-1));
	    cells.emplace_back(ij2kCell(i_f, (k_f-1)/2));
	  }
	}
	else // upper bound faces
	{
	  cells.emplace_back(ij2kCell(i_f-1, k_f));
	}
#endif

#ifdef __VALIDATION__
	if (!areVectorEqual(cells, cellsOfFace))
	{
	  cout << "getCellsOfFace: Nabla and index versions do not match !!\n";
	  cout << "  faceId = " << faceId << "\n";
	  cout << "  nabla = ";
	  for(int i(0) ; i < cellsOfFace.size() ; i++)
	    cout << cellsOfFace[i] << ", ";
	  cout << "\n";
	  cout << "  index = ";
	  for(int i(0) ; i < cells.size() ; i++)
	    cout << cells[i] << ", ";
	  cout << "\n";
	}
#endif

#ifdef __CARTESIAN__
	return cells;
#else
	return cellsOfFace;
#endif
}

vector<int>
CartesianMesh2D::getNeighbourCells(const int& cellId) const
{
#if defined(__VALIDATION__) || !defined(__CARTESIAN__)
	std::vector<int> neighbours;
	const auto& nodes(getNodesOfCell(cellId));
	for (auto nodeId : nodes)
	{
		auto adjacentCells(getCellsOfNode(nodeId));
		for (int quadId : adjacentCells)
			if (quadId != cellId)
				if (getNbCommonIds(nodes, m_geometry->getQuads()[quadId].getNodeIds()) == 2)
					neighbours.emplace_back(quadId);
	}
	std::sort(neighbours.begin(), neighbours.end());
	neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
#endif

	
#if defined(__VALIDATION__) || defined(__CARTESIAN__)
  	std::vector<int> neighbours_with_indices;
	int i,j;
	k2ijCell(cellId, i,j);
	if (i>=1) neighbours_with_indices.emplace_back(ij2kCell(i-1,j));
	if (i<m_nb_y_quads-1) neighbours_with_indices.emplace_back(ij2kCell(i+1,j));
	if (j>=1) neighbours_with_indices.emplace_back(ij2kCell(i,j-1));
	if (j<m_nb_x_quads-1) neighbours_with_indices.emplace_back(ij2kCell(i,j+1));
#endif
	  
#ifdef __VALIDATION__
	if (!areVectorEqual(neighbours, neighbours_with_indices))
	{
	  cout << "getCellsOfNode : Nabla and index versions do not match !!\n";
	  std::terminate();
	}
#endif

#ifdef __CARTESIAN__
	return neighbours_with_indices;
#else
	return neighbours;
#endif

}

vector<int>
CartesianMesh2D::getFacesOfCell(const int& cellId) const
{
#if defined(__VALIDATION__) || !defined(__CARTESIAN__)
	vector<int> cellEdgeIds;
	const auto& edges(m_geometry->getEdges());
	for (int edgeId=0; edgeId < edges.size(); ++edgeId)
		if (getNbCommonIds(edges[edgeId].getNodeIds(), m_geometry->getQuads()[cellId].getNodeIds()) == 2)
			cellEdgeIds.emplace_back(edgeId);
#endif

#if defined(__VALIDATION__) || defined(__CARTESIAN__)
	int i,j;
	k2ijCell(cellId, i,j);
	const int& bottom_face = 2*j + i * (2*m_nb_x_quads+1);
	const int& left_face = bottom_face + 1;
	const int& right_face = bottom_face + (j==m_nb_x_quads-1 ? 2 : 3);
	const int& top_face = bottom_face +  (i < m_nb_y_quads-1 ? 2*m_nb_x_quads+1 : 2*m_nb_x_quads+1 - j);
	vector<int> faces {bottom_face, left_face, right_face, top_face};
#endif

#ifdef __VALIDATION__
	if (!areVectorEqual(faces, cellEdgeIds))
	{
	  cout << "getFacesOfCell: Nabla and index versions do not match !!\n";
	  cout << "  cellId = " << cellId << "\n";
	  cout << "  nabla = ";
	  for(int i(0) ; i < cellEdgeIds.size() ; i++)
	    cout << cellEdgeIds[i] << ", ";
	  cout << "\n";
	  cout << "  index = ";
	  for(int i(0) ; i < faces.size() ; i++)
	    cout << faces[i] << ", ";
	  cout << "\n";
	  
	    

	  //std::terminate();
	}
#endif
	
#ifdef __CARTESIAN__
	return faces;
#else
	return cellEdgeIds;
#endif	
	
}

int
CartesianMesh2D::getCommonFace(const int& cellId1, const int& cellId2) const noexcept
{
	auto cell1Faces{getFacesOfCell(cellId1)};
	auto cell2Faces{getFacesOfCell(cellId2)};
	auto result = find_first_of(cell1Faces.begin(), cell1Faces.end(), cell2Faces.begin(), cell2Faces.end());
	if (result == cell1Faces.end())
	  return -1;
	else
	  return *result;
}

int
CartesianMesh2D::getBackCell(const int& faceId) const noexcept
{
  vector<int> cells = getCellsOfFace(faceId);
  if (cells.size() >= 1)
    return cells[0];
  else
    return -1;
}

int
CartesianMesh2D::getFrontCell(const int& faceId) const noexcept
{
  vector<int> cells = getCellsOfFace(faceId);
  if (cells.size() >= 2)
    return cells[1];
  else
    return -1;
}

  
int
CartesianMesh2D::getNbCommonIds(const vector<int>& as, const vector<int>& bs) const noexcept
{
	int nbCommonIds = 0;
	for (auto a : as)
		for (auto b : bs)
			if (a == b) nbCommonIds++;
	return nbCommonIds;
}

bool
CartesianMesh2D::isInnerEdge(const Edge& e) const noexcept
{
	return (find(m_inner_nodes.begin(), m_inner_nodes.end(), e.getNodeIds()[0]) != m_inner_nodes.end()) ||
	       (find(m_inner_nodes.begin(), m_inner_nodes.end(), e.getNodeIds()[1]) != m_inner_nodes.end());
}

bool
CartesianMesh2D::isInnerVerticalEdge(const Edge& e) const noexcept
{
	if (!isInnerEdge(e)) return false;
	return (e.getNodeIds()[0] == e.getNodeIds()[1] + m_nb_x_quads + 1 ||
			e.getNodeIds()[1] == e.getNodeIds()[0] + m_nb_x_quads + 1);
}

bool
CartesianMesh2D::isInnerHorizontalEdge(const Edge& e) const noexcept
{
	if (!isInnerEdge(e)) return false;
	return (e.getNodeIds()[0] == e.getNodeIds()[1] + 1 ||
			e.getNodeIds()[1] == e.getNodeIds()[0] + 1);
}


int
CartesianMesh2D::getBottomFaceOfCell(const int& cellId) const noexcept
{
  int i,j;
  k2ijCell(cellId, i,j);
  const int& bottom_face = 2*j + i * (2*m_nb_x_quads+1);
  return bottom_face;
}

int
CartesianMesh2D::getLeftFaceOfCell(const int& cellId) const noexcept
{
  int i,j;
  k2ijCell(cellId, i,j);
  const int& bottom_face = 2*j + i * (2*m_nb_x_quads+1);
  const int& left_face = bottom_face + 1;
  return left_face;
}

int
CartesianMesh2D::getRightFaceOfCell(const int& cellId) const noexcept
{
  int i,j;
  k2ijCell(cellId, i,j);
  const int& bottom_face = 2*j + i * (2*m_nb_x_quads+1);
  const int& right_face = bottom_face + (j==m_nb_x_quads-1 ? 2 : 3);
  return right_face;
}

int CartesianMesh2D::getTopFaceOfCell(const int& cellId) const noexcept
{
  int i,j;
  k2ijCell(cellId, i,j);
  const int& bottom_face = 2*j + i * (2*m_nb_x_quads+1);
  const int& top_face = bottom_face +  (i < m_nb_y_quads-1 ? 2*m_nb_x_quads+1 : 2*m_nb_x_quads+1 - j);
  return top_face;
}




}
