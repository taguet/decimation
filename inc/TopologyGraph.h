#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <set>
#include <map>
#include <limits>
#include "MeshUtils.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;
using Plane = Equation::Plane;
using Line = Equation::Line;

class TopologyGraph {
private:
	class Node {
	public:
		const int id;
		Plane plane;

		Node(TopologyGraph& parent, int id) : id{ id }, parent{ &parent } {}
		void add(Mesh::FaceHandle fh);
		void regroupIntoSelf(Node& region);
		const std::set<Mesh::FaceHandle>& getFaceHandles() const { return faces; }
		const std::set<Mesh::VertexHandle>& getVertexHandles() const { return vertices; }
		bool empty() const { return faces.empty(); }
		bool contains(Mesh::FaceHandle& fh) const { return faces.contains(fh); }
		float computeArea() const;
		void fitPlane();
		float sumVertexProjectedDistances() const;
	private:
		TopologyGraph* parent{ nullptr };
		std::set<Mesh::FaceHandle> faces{};
		std::set<Mesh::VertexHandle> vertices{};
	};

	std::map<int, Node> regions{};
	std::map<int, std::set<int>> edges{};
	Mesh& mesh;

	OpenMesh::FPropHandleT<int> f_group;

	int findTargetRegion(int regionID, float fitting_threshold) const;
	void regroupRegionIntoTarget(int regionID, int targetID);
	void ungroupRegion(int regionID, bool removeGroup);
	std::map<std::pair<int, int>, Line> filterLinesByRegion(const int region_id, const std::map<std::pair<int, int>, Line> &contour_lines) const;
	const Line& findClosestLine(const Mesh::EdgeHandle eh, const std::map<std::pair<int, int>, Line>& contour_lines) const;

public:
	const float area_threshold;
	const float fitting_threshold;

	TopologyGraph(Mesh& mesh, float area_threshold, float fitting_threshold);
	int size() {
		return regions.size();
	}

	void addFaceToRegion(int regionID, Mesh::FaceHandle fh);

	Node& getRegion(int regionID) { return regions.at(regionID); }
	const Node& getRegion(int regionID) const { return regions.at(regionID); }
	int getFaceRegion(Mesh::FaceHandle fh) const;
	Plane& getPlane(int regionID) { return getRegion(regionID).plane; }
	const Plane& getPlane(int regionID) const { return getRegion(regionID).plane; }
	std::set<int> getNeighbors(int regionID) const { return edges.at(regionID); }
	std::set<int> getRegionIDs() const;
	std::pair<int, int> getNeighborIDs(Mesh::EdgeHandle eh) const;

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
	bool areFacesInSameRegion(Mesh::FaceHandle fh_1, Mesh::FaceHandle fh_2) const;
	void fitPlanes();

	/// @brief Finds all intersections between each neighboring regions' planes
	/// @return The equation of a line.
	std::vector<Line> findPlanePlaneIntersections() const;
	/// @brief Find all intersections between region planes.
	/// @return An association between pairs of regions and lines.
	std::map<std::pair<int, int>, Line> findContourLines() const;
	/// @brief Finds all pairs of neighboring planes
	/// @return A container of region pairs.
	std::set<std::pair<const Plane*, const Plane*>> getNeighborPairs() const;
	/// @brief Finds all pairs of neighboring regions.
	/// @return A set of node pairs.
	std::set<std::pair<int, int>> getRegionPairs() const;

	/// @brief Finds the edges that belong to two different regions.
	/// @return A set of edges.
	std::set<Mesh::EdgeHandle> extractContour();

	void projectContourVertices(bool handleUngroupedFaces = true);

	int& faceGroup(Mesh::FaceHandle fh) {
		return mesh.property(f_group, fh);
	}
	int& faceGroup(Mesh::HalfedgeHandle heh) {
		return mesh.property(f_group, mesh.face_handle(heh));
	}
	int faceGroup(Mesh::FaceHandle fh) const {
		return mesh.property(f_group, fh);
	}
	int faceGroup(Mesh::HalfedgeHandle heh) const {
		return mesh.property(f_group, mesh.face_handle(heh));
	}
};