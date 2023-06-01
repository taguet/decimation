#include "TopologyGraph.h"

TopologyGraph::TopologyGraph(Mesh& mesh, float area_threshold, float fitting_threshold)
	: mesh{ mesh }, area_threshold{ area_threshold }, fitting_threshold{ fitting_threshold }
{}


void TopologyGraph::insertEdge(int node_1, int node_2) {
	if (node_1 == node_2)	return;
	std::set<int>& connected_nodes{ edges[node_1] };
	connected_nodes.insert(node_2);
}


void TopologyGraph::connectRegions(int regionID_1, int regionID_2) {
	if (regionID_1 == regionID_2)	return;
	insertEdge(regionID_1, regionID_2);
	insertEdge(regionID_2, regionID_1);
}


float TopologyGraph::Node::computeArea() const {
	float area{ 0.0f };
	for (auto face : faces) {
		area += MeshUtils::computeFaceArea(parent->mesh, face);
	}
	return area;
}


bool TopologyGraph::simplifyGraph() {
	for (auto it{ regions.begin() }; it != regions.end(); ) {
		Node& region{ it->second };
		if (region.computeArea() > area_threshold) {	//Check whether region is very small
			++it;
			continue;
		}
		else {
			int targetID{ findTargetRegion(region.id, 100.0f) }; //arbitrary threshold
			if (targetID == -1) {
				std::cout << "Deleted region " << region.id << std::endl;
				ungroupRegion(region.id);
				return true;
			}
			else {
				std::cout << "Merged region " << region.id << " into " << targetID << std::endl;
				regroupRegionIntoTarget(region.id, targetID);
				return true;
			}
		}
	}
	return false;
}


int TopologyGraph::findTargetRegion(int regionID, float fitting_threshold) {
	std::set<int> neighborIDs{ edges.at(regionID) };
	int targetID{ -1 };
	float min_sum{ std::numeric_limits<float>::max() };
	for (int id : neighborIDs) {
		Node& neighbor{ getRegion(id) };
		float sum{ neighbor.sumVertexProjectedDistances() };
		if (sum < min_sum) {
			min_sum = sum;
			targetID = id;
		}
	}
	if (min_sum <= fitting_threshold)
		return targetID;
	else
		return -1;
}


void TopologyGraph::fitPlanes() {
	for (auto& region_pair : regions) {
		region_pair.second.fitPlane();
	}
}


void TopologyGraph::regroupRegionIntoTarget(int regionID, int targetID) {
	Node& region{ getRegion(regionID) };
	Node& target{ getRegion(targetID) };
	target.regroupIntoSelf(region);
	std::set<int> neighborIDs{ edges.at(regionID) };
	for (int id : neighborIDs) {
		connectRegions(targetID, id);	//Target region inherits deleted region's edges
	}
	ungroupRegion(regionID);
}


void TopologyGraph::ungroupRegion(int regionID) {
	//Erase connections to removed region
	removeEdges(regionID);
	edges.erase(regionID);
	//Erase region node
	regions.erase(regionID);
}

void TopologyGraph::removeEdges(int regionID) {
	std::set<int> neighborIDs{ edges.at(regionID) };
	for (int id : neighborIDs) {
		edges.at(id).erase(regionID);
		edges.at(regionID).erase(id);
	}
}


void TopologyGraph::addFaceToRegion(int regionID, Mesh::FaceHandle fh) {
	auto it_region{ regions.find(regionID) };
	if (it_region == regions.end()) {
		regions.insert(std::make_pair(regionID, Node{ *this, regionID }));
	}
	regions.at(regionID).add(fh);
}


int TopologyGraph::getFaceRegion(Mesh::FaceHandle fh) {
	for (auto& p : regions) {
		Node& region{ p.second };
		if (region.contains(fh))
			return region.id;
	}
	return -1;
}


std::set<int> TopologyGraph::getRegionIDs() {
	std::set<int> ids{};
	for (auto& p : regions) {
		ids.insert(p.first);
	}
	return ids;
}


void TopologyGraph::Node::add(Mesh::FaceHandle fh) {
	Mesh& mesh{ parent->mesh };
	faces.insert(fh);
	Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(fh) };
	vertices.insert(mesh.from_vertex_handle(heh));
	vertices.insert(mesh.to_vertex_handle(heh));
	vertices.insert(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
}


void TopologyGraph::Node::regroupIntoSelf(Node& region) {
	faces.insert(region.faces.begin(), region.faces.end());
	vertices.insert(region.vertices.begin(), region.vertices.end());
	fitPlane();
}


void TopologyGraph::Node::fitPlane() {
	std::set<Mesh::VertexHandle> vertices{};
	plane_params = MeshUtils::fitPlaneToVertices(parent->mesh, vertices);
}


float TopologyGraph::Node::sumVertexProjectedDistances() {
	Mesh& mesh{ parent->mesh };
	float sum{ 0.0f };
	for (auto& vh : vertices) {
		Mesh::Point p{ mesh.point(vh) };
		float a{ plane_params[1] };
		float b{ plane_params[2] };
		float c{ -1.0f };	// Ax + By + D = Z
		float d{ plane_params[0] };
		sum += std::abs(a * p[0] + b * p[1] + c * p[2] + d) / std::sqrtf(a * a + b * b + c * c);
	}
	return sum;
}
