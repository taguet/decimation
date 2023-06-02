#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <set>
#include <map>
#include <limits>
#include "MeshUtils.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class TopologyGraph {
private:
	class Node {
	public:
		const int id;
		Vector3f plane_params;

		Node(TopologyGraph& parent, int id) : id{ id }, parent{ &parent } {}
		void add(Mesh::FaceHandle fh);
		void regroupIntoSelf(Node& region);
		const std::set<Mesh::FaceHandle>& getFaceHandles() { return faces; }
		const std::set<Mesh::VertexHandle>& getVertexHandles() { return vertices; }
		bool empty() const { return faces.empty(); }
		bool contains(Mesh::FaceHandle& fh) { return faces.find(fh) != faces.end(); }
		float computeArea() const;
		void fitPlane();
		float sumVertexProjectedDistances();
	private:
		TopologyGraph* parent{ nullptr };
		std::set<Mesh::FaceHandle> faces{};
		std::set<Mesh::VertexHandle> vertices{};
	};

	std::map<int, Node> regions{};
	std::map<int, std::set<int>> edges{};
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
	int getFaceRegion(Mesh::FaceHandle fh);
	std::set<int> getNeighbors(int regionID) { return edges.at(regionID); }
	std::set<int> getRegionIDs();

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