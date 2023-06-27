#include "MeshTools.h"
#include <list>

void MeshTools::extractRegions(TopologyGraph& graph) {
	//initialization
	std::list<OpenMesh::FaceHandle> ungrouped_faces{};
	for (auto f_iter{ mesh_->faces_begin() }; f_iter != mesh_->faces_end(); ++f_iter) {
		graph.faceGroup(f_iter) = -1;
		ungrouped_faces.push_back(f_iter);
	}
	std::cerr << "Executing region growing algorithm...\n";
	growRegions(ungrouped_faces, graph);
	std::cerr << "Done\nBuilding topology graph...\n";
	buildTopologyGraph(graph);
}


void MeshTools::growRegions(std::list<Mesh::FaceHandle>& ungrouped_faces, TopologyGraph& graph) {
	AnisotropicLaplacian lapl{ *mesh_ };	//We need this for filtered normals
	int g{-1};
	while (!ungrouped_faces.empty()) {
		Mesh::FaceHandle fh{ ungrouped_faces.front() };
		ungrouped_faces.pop_front();
		graph.faceGroup(fh) = ++g;
		graph.addFaceToRegion(g, fh);
		std::list<Mesh::FaceHandle>  neighbors{  };
		MeshUtils::getFaceNeighbors(*mesh_, fh, neighbors);
		std::cerr << "Region " << g << '\n';
		for (auto f_neighbor{ neighbors.begin() }; f_neighbor != neighbors.end(); ) {
			Mesh::Normal f_normal_neighbor{ lapl.filterFaceNormal(*f_neighbor) };
			Mesh::Normal f_normal{ lapl.filterFaceNormal(fh) };
			if (graph.faceGroup(*f_neighbor) != -1) {
				f_neighbor = neighbors.erase(f_neighbor);
				continue;
			}
			else if (dot(f_normal_neighbor, f_normal) > sin(0.349066f)) {
				graph.faceGroup(*f_neighbor) = g;
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
			if (graph.faceGroup(f_iter) != graph.faceGroup(neighbor)) {
				graph.connectRegions(graph.faceGroup(f_iter), graph.faceGroup(neighbor));
			}
		}
	}
	std::cerr << "Done\nFitting planes...\n";
	graph.fitPlanes();
	std::cerr << "Done\nSimplifying graph...\n";
	while (graph.simplifyGraph());	//simplify graph until no changes are made
	std::cerr << "Done\n";
}
