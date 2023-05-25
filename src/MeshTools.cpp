#include "MeshTools.h"
#include <list>
#include <limits>

void MeshTools::extractRegions() {
	//initialization
	mesh_->add_property(f_group);
	std::list<OpenMesh::FaceHandle> ungrouped_faces{};
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		faceGroup(f_iter) = -1;
		ungrouped_faces.push_back(f_iter);
	}

	TopologyGraph graph{ *this, 1.0f, 1.0f};
	growRegions(ungrouped_faces, graph);
	buildTopologyGraph(graph);
	//TODO delete following lines
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		faceGroup(f_iter) = -1;
	}
	graph.updateIndices();
}


void MeshTools::TopologyGraph::updateIndices() {
	for (auto& id_region_pair : regions) {
		id_region_pair.second.updateIndices();
	}
}

void MeshTools::TopologyGraph::Node::updateIndices() {
	for (auto& face : faces) {
		parent->parent->faceGroup(face) = id;
	}
}


void MeshTools::growRegions(std::list<Mesh::FaceHandle>& ungrouped_faces, TopologyGraph& graph) {
	AnisotropicLaplacian lapl{ *mesh_ };	//We need this for filtered normals
	int g{-1};
	while (!ungrouped_faces.empty()) {
		Mesh::FaceHandle fh{ ungrouped_faces.front() };
		ungrouped_faces.pop_front();
		faceGroup(fh) = ++g;
		graph.addFaceToRegion(g, fh);
		std::list<Mesh::FaceHandle>  neighbors{  };
		MeshUtils::getFaceNeighbors(*mesh_, fh, neighbors);
		for (auto f_neighbor{ neighbors.begin() }; f_neighbor != neighbors.end(); ) {
			Mesh::Normal f_normal_neighbor{ lapl.filterFaceNormal(*f_neighbor) };
			Mesh::Normal f_normal{ lapl.filterFaceNormal(fh) };
			if (faceGroup(*f_neighbor) != -1) {
				f_neighbor = neighbors.erase(f_neighbor);
				continue;
			}
			else if (dot(f_normal_neighbor, f_normal) > sin(0.349066f)) {
				faceGroup(*f_neighbor) = g;
				graph.addFaceToRegion(g, *f_neighbor);
				std::list<Mesh::FaceHandle> extended_neighborhood{  };
				MeshUtils::getFaceNeighbors(*mesh_, *f_neighbor, extended_neighborhood);
				ungrouped_faces.remove(*f_neighbor);
				neighbors.splice(neighbors.end(), extended_neighborhood);
			}
			++f_neighbor;
		}
	}
}


void MeshTools::buildTopologyGraph(TopologyGraph& graph) {
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		for (auto neighbor{ mesh_->ff_iter(f_iter) }; neighbor; ++neighbor) {
			if (faceGroup(f_iter) != faceGroup(neighbor)) {
				graph.connectRegions(faceGroup(f_iter), faceGroup(neighbor));
			}
		}
	}
	graph.fitPlanes();
	while (graph.simplifyGraph());	//simplify graph until no changes are made
}


Vector3f MeshTools::fitPlaneToVertices(std::set<Mesh::VertexHandle>& vertices) const {
	MatrixXf regressors{ vertices.size(), 3 };	// x & y
	VectorXf observed{ vertices.size() };	// z
	Vector3f parameters{ 0.0f, 0.0f, 0.0f };	//c, a & b

	{
		std::set<Mesh::VertexHandle>::iterator it{ vertices.begin() };
		for (int i{ 0 };	i < regressors.rows(), it != vertices.end(); ++i, ++it) {
			Mesh::Point p{ mesh_->point(*it) };
			regressors(i, 0) = 1.0f;
			regressors(i, 1) = p[0];
			regressors(i, 2) = p[1];
			observed(i) = p[2];
		}
	}
	MatrixXf t_regressors{ regressors.transpose() };
	parameters = (t_regressors * regressors).inverse() * t_regressors * observed;
	return parameters;
}


MeshTools::TopologyGraph::TopologyGraph(MeshTools& parent, float area_threshold, float fitting_threshold) 
	: parent{ &parent}, area_threshold{ area_threshold }, fitting_threshold{fitting_threshold}
{}


void MeshTools::TopologyGraph::insertEdge(int node_1, int node_2) {
	if (node_1 == node_2)	return;
	std::set<int>& connected_nodes{ edges[node_1] };
	connected_nodes.insert(node_2);
}


void MeshTools::TopologyGraph::connectRegions(int regionID_1, int regionID_2) {
	if (regionID_1 == regionID_2)	return;
	insertEdge(regionID_1, regionID_2);
	insertEdge(regionID_2, regionID_1);
}


float MeshTools::TopologyGraph::Node::computeArea() const {
	float area{ 0.0f };
	for (auto face : faces) {
		area += MeshUtils::computeFaceArea(*parent->parent->mesh_, face);
	}
	return area;
}


bool MeshTools::TopologyGraph::simplifyGraph() {
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


int MeshTools::TopologyGraph::findTargetRegion(int regionID, float fitting_threshold)  {
	std::set<int> neighborIDs{ edges.at(regionID)};
	int targetID{ -1 };
	float min_sum{ numeric_limits<float>::max() };
	for (int id : neighborIDs) {
		Node& neighbor{ getRegion(id)};
		float sum{ neighbor.sumVertexProjectedDistances()};
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


void MeshTools::TopologyGraph::fitPlanes() {
	for (auto& region_pair : regions) {
		region_pair.second.fitPlane();
	}
}


void MeshTools::TopologyGraph::regroupRegionIntoTarget(int regionID, int targetID) {
	Node& region{ getRegion(regionID) };
	Node& target{ getRegion(targetID) };
	target.regroupIntoSelf(region);
	std::set<int> neighborIDs{ edges.at(regionID) };
	for (int id : neighborIDs) {
		connectRegions(targetID, id);	//Target region inherits deleted region's edges
	}
	ungroupRegion(regionID);
}


void MeshTools::TopologyGraph::ungroupRegion(int regionID) {
	//Erase connections to removed region
	removeEdges(regionID);
	edges.erase(regionID);
	//Erase region node
	regions.erase(regionID);
}

void MeshTools::TopologyGraph::removeEdges(int regionID) {
	std::set<int> neighborIDs{ edges.at(regionID) };
	for (int id : neighborIDs) {
		edges.at(id).erase(regionID);
		edges.at(regionID).erase(id);
	}
}


void MeshTools::TopologyGraph::addFaceToRegion(int regionID, Mesh::FaceHandle fh) {
	auto it_region{ regions.find(regionID) };
	if (it_region == regions.end()) {
		regions.insert(std::make_pair(regionID, Node{ *this, regionID }));
	}
	regions.at(regionID).add(fh);
}


void MeshTools::TopologyGraph::Node::add(Mesh::FaceHandle fh) {
	Mesh* mesh{ parent->parent->mesh_ };
	faces.insert(fh); 
	Mesh::HalfedgeHandle heh{ mesh->halfedge_handle(fh) };
	vertices.insert(mesh->from_vertex_handle(heh));
	vertices.insert(mesh->to_vertex_handle(heh));
	vertices.insert(mesh->to_vertex_handle(mesh->next_halfedge_handle(heh)));
}


void MeshTools::TopologyGraph::Node::regroupIntoSelf(Node& region) {
	faces.insert(region.faces.begin(), region.faces.end());
	vertices.insert(region.vertices.begin(), region.vertices.end());
	fitPlane();
}


void MeshTools::TopologyGraph::Node::fitPlane() {
	std::set<Mesh::VertexHandle> vertices{};
	plane_params = parent->parent->fitPlaneToVertices(vertices);
}


float MeshTools::TopologyGraph::Node::sumVertexProjectedDistances() {
	Mesh* mesh{ parent->parent->mesh_ };
	float sum{ 0.0f };
	for (auto& vh : vertices) {
		Mesh::Point p{ mesh->point(vh) };
		float a{ plane_params[1] };
		float b{ plane_params[2] };
		float c{ 1.0f };
		float d{ plane_params[0] };
		sum += std::abs(a * p[0] + b * p[1] + c * p[2] + d) / std::sqrtf(a * a + b * b + c * c);
	}
	return sum;
}