#include "TopologyGraph.h"

TopologyGraph::TopologyGraph(Mesh& mesh, float area_threshold, float fitting_threshold)
	: mesh{ mesh }, area_threshold{ area_threshold }, fitting_threshold{ fitting_threshold }
{
	this->mesh.add_property(f_group);
	for (auto f_it{ mesh.faces_begin() }; f_it != mesh.faces_end(); ++f_it) {
		faceGroup(f_it) = -1;
	}
}


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
				ungroupRegion(region.id, true);
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
		connectRegions(targetID, id);	//Target region inherits deleted region's neighbors
	}
	ungroupRegion(regionID, false);
}


void TopologyGraph::ungroupRegion(int regionID, bool removeGroup) {
	if (removeGroup) {
		Node& region{ regions.find(regionID)->second };
		for (auto fh : region.getFaceHandles()) {
			faceGroup(fh) = -1;
		}
	}
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


std::set<int> TopologyGraph::getRegionIDs() const {
	std::set<int> ids{};
	for (auto& p : regions) {
		ids.insert(p.first);
	}
	return ids;
}


bool TopologyGraph::areFacesInSameRegion(Mesh::FaceHandle fh_1, Mesh::FaceHandle fh_2) {
	return faceGroup(fh_1) == faceGroup(fh_2);
}


std::set<Mesh::EdgeHandle> TopologyGraph::extractContour() {
	std::set<Mesh::EdgeHandle> contour_edges{};
	for (auto f_it{ mesh.faces_begin() }; f_it != mesh.faces_end(); ++f_it) {
		for (auto fh_it{ mesh.fh_iter(f_it)}; fh_it; ++fh_it) {
			Mesh::FaceHandle neighbor_face{ mesh.face_handle(mesh.opposite_halfedge_handle(fh_it)) };
			if (!areFacesInSameRegion(f_it, neighbor_face)) {
				contour_edges.insert(mesh.edge_handle(fh_it));
			}
		}
	}
	return contour_edges;
}


std::vector<Line> TopologyGraph::findPlanePlaneIntersections() {
	std::set<std::pair<Plane*, Plane*>>& neighbor_pairs{ getNeighborPairs() };
	std::vector<Line> lines;
	lines.reserve(neighbor_pairs.size());
	for (auto const & [plane_1, plane_2] : neighbor_pairs) {
		lines.push_back(plane_1->findPlanePlaneIntersection(*plane_2));
	}
	return lines;
}


std::map<std::pair<int, int>, Line> TopologyGraph::findContourLines() {
	std::map<std::pair<int, int>, Line> lines;
	std::set<std::pair<int, int>> region_pairs{ getRegionPairs() };
	for (auto const & [id_1, id_2] : region_pairs) {
		if (!lines.contains({ id_1, id_2 })) {
			Plane* plane_1{ &getRegion(id_1).plane };
			Plane* plane_2{ &getRegion(id_2).plane };
			Line line{ plane_1->findPlanePlaneIntersection(*plane_2) };
			lines.insert({ {id_1, id_2}, line });
		}
	}
	return lines;
}


std::set<std::pair<int, int>> TopologyGraph::getRegionPairs() {
	std::set<std::pair<int, int>> region_pairs;
	for (auto const & [id, neighbors] : edges) {
		for (auto const& neighbor_id : neighbors) {
			std::pair<int, int> region_pair{ id, neighbor_id };
			std::pair<int, int> reverse{ neighbor_id, id };
			if (!region_pairs.contains(reverse)) {
				region_pairs.insert(region_pair);
			}
		}
	}
	return region_pairs;
}


std::set<std::pair<Plane*, Plane*>> TopologyGraph::getNeighborPairs() {
	std::set<std::pair<Plane*, Plane*>> plane_pairs{};
	for (auto const & [id, neighbors] : edges) {
		for (int neighbor_id : neighbors) {
			Node& region_1{ getRegion(id) };
			Node& region_2{ getRegion(neighbor_id) };
			std::pair<Plane*, Plane*> plane_pair{ &region_1.plane, &region_2.plane };
			std::pair<Plane*, Plane*> rev{ plane_pair.second, plane_pair.first };
			if (!plane_pairs.contains(rev)) {
				plane_pairs.insert(plane_pair);
			}
		}
	}
	return plane_pairs;
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
	for (auto fh : region.getFaceHandles()) {
		parent->faceGroup(fh) = id;
	}
	fitPlane();
}


void TopologyGraph::Node::fitPlane() {
	plane = MeshUtils::fitPlaneToVertices(parent->mesh, vertices);
}


float TopologyGraph::Node::sumVertexProjectedDistances() {
	Mesh& mesh{ parent->mesh };
	float sum{ 0.0f };
	for (auto& vh : vertices) {
		Mesh::Point p{ mesh.point(vh) };
		sum += plane.distToPoint(Vector3f{ p[0], p[1], p[2]});
	}
	return sum;
}
