#include "MeshTools.h"
#include <list>

void MeshTools::extractRegions() {
	//initialization
	mesh_->add_property(f_group);
	std::list<OpenMesh::FaceHandle> ungrouped_faces{};
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		faceGroup(f_iter) = -1;
		ungrouped_faces.push_back(f_iter);
	}

	TopologyGraph graph{*mesh_, 1.0f, 1.0f};
	growRegions(ungrouped_faces, graph);
	buildTopologyGraph(graph);
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
				graph.insertEdge(faceGroup(f_iter), faceGroup(neighbor));
				graph.insertEdge(faceGroup(neighbor), faceGroup(f_iter));
			}
		}
	}
}


MeshTools::TopologyGraph::TopologyGraph(Mesh& mesh_, float area_threshold, float fitting_threshold) 
	: mesh_{ &mesh_ }, area_threshold{ area_threshold }, fitting_threshold{fitting_threshold}
{}


void MeshTools::TopologyGraph::insertEdge(int node_1, int node_2) {
	std::set<int>& connected_nodes{ edges[node_1] };
	connected_nodes.insert(node_2);
}


float MeshTools::TopologyGraph::Node::computeArea() const {
	float area{ 0.0f };
	for (auto face : faces) {
		area += MeshUtils::computeFaceArea(*parent->mesh_, face);
	}
	return area;
}


void MeshTools::TopologyGraph::simplifyGraph() {
	for (auto it{ regions.begin() }; it != regions.end(); ) {
		Node& region{ it->second };
		if (region.computeArea() > area_threshold) {	//Check whether region is very small
			++it;
			continue;
		}
		else {
			//TODO
		}
	}
}


int MeshTools::TopologyGraph::findTargetRegion(int regionID) const {
	return -1; //TODO
}


VectorXd MeshTools::fitPlaneToVertices(std::set<Mesh::VertexHandle>& vertices) const {
	MatrixXd regressors{ vertices.size(), 3 };
	VectorXd observed{ vertices.size() };
	VectorXd parameters{ vertices.size() };

	for (std::set<Mesh::VertexHandle>::iterator it{vertices.begin()}, int i{0}; 
		i < regressors.rows(), it != vertices.end(); ++i, ++it) {
		Mesh::Point p{ mesh_->point(*it)};
		regressors(i, 0) = 1;
		regressors(i, 1) = p[0];	//x
		regressors(i, 2) = p[1];	//y
		observed(i) = p[2];			//z
	}
	MatrixXd t_regressors{ regressors.transpose() };
	parameters = (t_regressors * regressors).inverse() * t_regressors * observed;
	return parameters;
}