#pragma once

#include "gl.h"
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "Laplacian.h"


typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class MeshTools
{
private:
	Mesh* mesh_{ nullptr };
public:
	
	MeshTools(void);
	MeshTools(Mesh& mesh);
	void setMesh(Mesh &mesh);


	template <typename T>
	void taubinSmoothing(int iterations = 1, float lambda = 0.2, float mu = -0.2) {
		std::unique_ptr<Laplacian> laplacian{ new T(*mesh_) };
		for (int i{ 0 }; i < iterations; ++i) {
			laplacian->computeLaplacians();
			for (auto v_it{ mesh_->vertices_begin() }; v_it != mesh_->vertices_end(); ++v_it) {
				mesh_->set_point(v_it, mesh_->point(v_it) + laplacian->laplacian_displacement(v_it) * lambda);
			}
			mesh_->update_normals();

			laplacian->computeLaplacians();
			for (auto v_it{ mesh_->vertices_begin() }; v_it != mesh_->vertices_end(); ++v_it) {
				mesh_->set_point(v_it, mesh_->point(v_it) + laplacian->laplacian_displacement(v_it) * mu);
			}
			mesh_->update_normals();
		}
	}

	/// Smoothes mesh using Laplacian smoothing. Template type should derive Laplacian class.
	/// iterations is the number of iteration of the given laplacian algorithm
	/// factor is the scalar diffusion coefficient
	template <typename T>
	void smoothMesh(int iterations=1, float factor=0.2) {
		std::unique_ptr<Laplacian> laplacian{ new T(*mesh_) };
		for (int i{ 0 }; i < iterations; ++i) {
			std::cout << "Iteration " << i + 1 << std::endl;
			laplacian->computeLaplacians();
			for (auto v_it{ mesh_->vertices_begin() }; v_it != mesh_->vertices_end(); ++v_it) {
				mesh_->set_point(v_it, mesh_->point(v_it) + laplacian->laplacian_displacement(v_it) * factor);
			}
			mesh_->update_normals();
		}
	}
};

