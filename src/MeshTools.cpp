#include "MeshTools.h"
#include <math.h>


MeshTools::MeshTools(void) {
	
}

MeshTools::MeshTools(Mesh& mesh) {
	this->mesh_ = &mesh;
}

void MeshTools::setMesh(Mesh& mesh) {
	this->mesh_ = &mesh;
}


/*
void MeshTools::taubinSmoothing(float lambda, float mu, int iterations) {
	for (int i{ 0 }; i < iterations; ++i) {
		variance_edge_length = computeVarianceEdgeLength();
		calc_discrete_laplacian();
		for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
			mesh_.set_point(v_it, mesh_.point(v_it) + laplacian_displacement(v_it) * lambda);
		}
		variance_edge_length = computeVarianceEdgeLength();
		mesh_.update_normals();
		calc_discrete_laplacian();
		for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
			mesh_.set_point(v_it, mesh_.point(v_it) + laplacian_displacement(v_it) * mu);
		}
		mesh_.update_normals();
	}
}
*/
void MeshTools::smoothMesh(Laplacian& laplacian, int iterations, float factor) {
	for (int i{ 0 }; i < iterations; ++i) {
		std::cout << "Iteration " << i + 1 << std::endl;
		laplacian.computeLaplacians();
		for (auto v_it{ mesh_->vertices_begin() }; v_it != mesh_->vertices_end(); ++v_it) {
			//std::cout << "Point: " << mesh_.point(v_it) << "\tLaplacien: " << laplacian_displacement(v_it) << std::endl;
			mesh_->set_point(v_it, mesh_->point(v_it) + laplacian.laplacian_displacement(v_it) * factor);
		}
		mesh_->update_normals();
	}
}