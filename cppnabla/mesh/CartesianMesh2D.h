/*******************************************************************************
 * Copyright (c) 2020 CEA
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 * Contributors: see AUTHORS file
 *******************************************************************************/
#ifndef MESH_CARTESIANMESH2D_H_
#define MESH_CARTESIANMESH2D_H_

#include "types/Types.h"
#include "mesh/MeshGeometry.h"

using namespace std;

namespace nablalib
{

class CartesianMesh2D
{
public:
	static constexpr int MaxNbNodesOfCell = 4;
	static constexpr int MaxNbNodesOfFace = 2;
	static constexpr int MaxNbCellsOfNode = 4;
	static constexpr int MaxNbCellsOfFace = 2;
	static constexpr int MaxNbFacesOfCell = 4;
	static constexpr int MaxNbNeighbourCells = 2;

	CartesianMesh2D(MeshGeometry<2>* geometry,
		const vector<int>& inner_nodes_ids,
		const vector<int>& top_nodes_ids,
		const vector<int>& bottom_nodes_ids,
		const vector<int>& left_nodes_ids,
		const vector<int>& right_nodes_ids,
		const int top_left_node_id,
		const int top_right_node_id,
		const int bottom_left_node_id,
		const int bottom_right_node_id);

	MeshGeometry<2>* getGeometry() noexcept { return m_geometry; }

	size_t getNbNodes() const noexcept { return m_geometry->getNodes().size(); }
	size_t getNbCells() const noexcept { return m_geometry->getQuads().size(); }
	size_t getNbFaces() const noexcept { return m_geometry->getEdges().size(); }
	vector<int> getFaces() const noexcept { return m_faces; }
	
	size_t getNbInnerNodes() const noexcept { return m_inner_nodes.size(); }
	const vector<int>& getInnerNodes() const noexcept { return m_inner_nodes; }

	size_t getNbTopNodes() const noexcept { return m_top_nodes.size(); }
	const vector<int>& getTopNodes() const noexcept { return m_top_nodes; }

	size_t getNbBottomNodes() const noexcept { return m_bottom_nodes.size(); }
	const vector<int>& getBottomNodes() const noexcept { return m_bottom_nodes; }

	size_t getNbLeftNodes() const noexcept { return m_left_nodes.size(); }
	const vector<int>& getLeftNodes() const noexcept { return m_left_nodes; }

	size_t getNbRightNodes() const noexcept { return m_right_nodes.size(); }
	const vector<int>& getRightNodes() const noexcept { return m_right_nodes; }

	size_t getNbOuterFaces() const noexcept { return m_outer_faces.size(); }
	vector<int> getOuterFaces() const noexcept { return m_outer_faces; }

	size_t getNbInnerFaces() const noexcept { return m_inner_faces.size(); }
	vector<int> getInnerFaces() const noexcept { return m_inner_faces; }
	
	size_t getNbInnerHorizontalFaces() const noexcept { return m_inner_horizontal_faces.size(); }
	vector<int> getInnerHorizontalFaces() const noexcept { return m_inner_horizontal_faces; }

	size_t getNbInnerVerticalFaces() const noexcept { return m_inner_vertical_faces.size(); }
	vector<int> getInnerVerticalFaces() const noexcept { return m_inner_vertical_faces; }

	size_t getNbTopLeftNode() const noexcept { return 1; }
	vector<int> getTopLeftNode() const noexcept { return vector<int>({ m_top_left_node}); }

	size_t getNbTopRightNode() const noexcept { return 1; }
	vector<int> getTopRightNode() const noexcept { return vector<int>({ m_top_right_node}); }

	size_t getNbBottomLeftNode() const noexcept { return 1; }
	vector<int> getBottomLeftNode() const noexcept { return vector<int>({ m_bottom_left_node}); }

	size_t getNbBottomRightNode() const noexcept { return 1; }
	vector<int> getBottomRightNode() const noexcept { return vector<int>({ m_bottom_right_node}); }
		
	const array<int, 4>& getNodesOfCell(const int& cellId) const noexcept;
	const array<int, 2>& getNodesOfFace(const int& faceId) const noexcept;
	vector<int> getCellsOfNode(const int& nodeId) const noexcept;
	vector<int> getCellsOfFace(const int& faceId) const;
	vector<int> getNeighbourCells(const int& cellId) const;
	vector<int> getFacesOfCell(const int& cellId) const;
	int getCommonFace(const int& cellId1, const int& cellId2) const noexcept;
	int getBackCell(const int& faceId) const noexcept;
	int getFrontCell(const int& faceId) const noexcept;
	int getFirstNodeOfFace(const int& faceId) const noexcept;
	int getSecondNodeOfFace(const int& faceId) const noexcept;

	int getBottomFaceOfCell(const int& cellId) const noexcept;
	int getLeftFaceOfCell(const int& cellId) const noexcept;
	int getRightFaceOfCell(const int& cellId) const noexcept;
	int getTopFaceOfCell(const int& cellId) const noexcept;

	
	
	inline int  ij2kCell(const int& i, const int&j) const noexcept;
	inline void k2ijCell(const int& k, int& i, int&j) const noexcept;
	inline int  ij2kNode(const int& i, const int&j) const noexcept;
	inline void k2ijNode(const int& k, int& i, int&j) const noexcept;
	
private:
	MeshGeometry<2>* m_geometry;
	vector<int> m_inner_nodes;
	vector<int> m_top_nodes;
	vector<int> m_bottom_nodes;
	vector<int> m_left_nodes;
	vector<int> m_right_nodes;
	vector<int> m_outer_faces;
	vector<int> m_faces;
	vector<int> m_inner_horizontal_faces;
	vector<int> m_inner_vertical_faces;
	vector<int> m_inner_faces;
	int m_top_left_node;
	int m_top_right_node;
	int m_bottom_left_node;
	int m_bottom_right_node;

	int m_nb_x_quads;
	int m_nb_y_quads;

	bool isInnerEdge(const Edge& e) const noexcept;
	bool isInnerVerticalEdge(const Edge& e) const noexcept;
	bool isInnerHorizontalEdge(const Edge& e) const noexcept;

	int getNbCommonIds(const vector<int>& a, const vector<int>& b) const noexcept;
	template <std::size_t T, std::size_t U>
	int	getNbCommonIds(const std::array<int, T>& as, const std::array<int, U>& bs) const noexcept
	{
		int nbCommonIds(0);
		for (const auto& a : as)
			if (find(bs.begin(), bs.end(), a) != bs.end())
				++nbCommonIds;
		return nbCommonIds;
	}
};

}
#endif /* MESH_CARTESIANMESH2D_H_ */
