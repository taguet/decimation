#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <set>
#include <map>
#include <Eigen/Dense>
#include <limits>
#include "MeshUtils.h"

using Eigen::Vector3f;

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class TopologyGraph {
private:
	class Node {
	public:
		const int id;

		Node(TopologyGraph& parent, int id) : id{ id }, parent{ &parent } {}
		void add(Mesh::FaceHandle fh);
		void regroupIntoSelf(Node& region);
		bool empty() const { return faces.empty(); }
		float computeArea() const;
		void fitPlane();
		float sumVertexProjectedDistances();
	private:
		TopologyGraph* parent{ nullptr };
		std::set<Mesh::FaceHandle> faces;
		std::set<Mesh::VertexHandle> vertices;
		Vector3f plane_params;
	};

	std::map<int, Node> regions;
	std::map<int, std::set<int>> edges;
	Mesh& mesh;

	int findTargetRegion(int regionID, float fitting_threshold);
	void regroupRegionIntoTarget(int regionID, int targetID);
	void ungroupRegion(int regionID);

public:
	const float area_threshold;
	const float fitting_threshold;

	TopologyGraph(Mesh& mesh, float area_threshold, float fitting_threshold);
	int size() {
		return regions.size();
	}

	void addFaceToRegion(int regionID, Mesh::FaceHandle fh);

	Node& getRegion(int regionID) { return regions.at(regionID); }

	/// @brief Insert an edge from node_1 to node_2.
	/// @param node_1 Start node.
	/// @param node_2 End node.
	void insertEdge(int node_1, int node_2);

	/// @brief Connect two region nodes with an edge, going both ways.
	/// @param regionID_1 
	/// @param regionID_2 
	void connectRegions(int regionID_1, int regionID_2);

	void removeEdges(int regionID);

	bool simplifyGraph();
	void fitPlanes();
};