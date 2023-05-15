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

	growRegions(ungrouped_faces);
	//TODO build graph
}


void MeshTools::growRegions(std::list<Mesh::FaceHandle>& ungrouped_faces) {
	AnisotropicLaplacian lapl{ *mesh_ };	//We need this for filtered normals
	int g{-1};
	while (!ungrouped_faces.empty()) {
		Mesh::FaceHandle fh{ ungrouped_faces.front() };
		ungrouped_faces.pop_front();
		faceGroup(fh) = ++g;
		std::list<Mesh::FaceHandle>  neighbors{  };
		neighbors = MeshUtils::getFaceNeighbors(*mesh_, fh, neighbors);
		for (auto f_neighbor{ neighbors.begin() }; f_neighbor != neighbors.end(); ) {
			Mesh::Normal f_normal_neighbor{ lapl.filterFaceNormal(*f_neighbor) };
			Mesh::Normal f_normal{ lapl.filterFaceNormal(fh) };
			if (faceGroup(*f_neighbor) != -1) {
				f_neighbor = neighbors.erase(f_neighbor);
				continue;
			}
			else if (dot(f_normal_neighbor, f_normal) > sin(0.349066f)) {
				faceGroup(*f_neighbor) = g;
				std::list<Mesh::FaceHandle> extended_neighborhood{  };
				extended_neighborhood = MeshUtils::getFaceNeighbors(*mesh_, *f_neighbor, extended_neighborhood);
				ungrouped_faces.remove(*f_neighbor);
				neighbors.splice(neighbors.end(), extended_neighborhood);
			}
			++f_neighbor;
		}
	}
}