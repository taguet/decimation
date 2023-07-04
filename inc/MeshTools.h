#pragma once

#include "gl.h"
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <map>
#include <set>
#include "Laplacian.h"
#include "TopologyGraph.h""


typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

using Quadric = Eigen::Matrix4f;

/// @brief Toolset class for modifying a mesh.
class MeshTools
{
private:
	Mesh* mesh_{ nullptr };
public:
	
	MeshTools(void) = default;
	MeshTools(Mesh& mesh) { this->mesh_ = &mesh; }
	void setMesh(Mesh &mesh) { this->mesh_ = &mesh; }

	/// @brief Smoothes mesh with Taubin smoothing.
	/// 
	/// Avoids shrinkage of mesh by growing the mesh after each step.
	/// 
	/// @tparam T Derived from the Laplacian class.
	/// @param iterations 
	/// @param lambda Shrinking coefficient. Should be positive.
	/// @param mu Growing coefficient. Should be negative.
	template <typename T>
	void taubinSmoothing(int iterations = 1, float lambda = 0.2, float mu = -0.2) {
		static_assert(std::is_base_of_v<Laplacian, T>);
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

	/// @brief Smoothes mesh using Laplacian smoothing.
	/// @tparam T Derived from the Laplacian class.
	/// @param iterations Number or iterations over the mesh.
	/// @param factor A scalar diffusion factor. Give it a lower value if the smoothing results in a spiky mesh.
	template <typename T>
	void smoothMesh(int iterations=1, float factor=0.2) {
		static_assert(std::is_base_of_v<Laplacian, T>);
		std::unique_ptr<Laplacian> laplacian{ new T(*mesh_) };
		for (int i{ 0 }; i < iterations; ++i) {
			std::cerr << "Iteration " << i + 1 << '\t';
			laplacian->computeLaplacians();
			for (auto v_it{ mesh_->vertices_begin() }; v_it != mesh_->vertices_end(); ++v_it) {
				mesh_->set_point(v_it, mesh_->point(v_it) + laplacian->laplacian_displacement(v_it) * factor);
			}
			mesh_->update_normals();
			std::cerr << "OK\n";
		}
	}


	void extractRegions(TopologyGraph& graph);

private:
	/// @brief Region growing algorithm to detect and separate all planar regions.
	/// @param ungrouped_faces Every faces in the mesh that haven't been sorted in a group yet.
	/// @param graph The topology graph to update throughout the process.
	void growRegions(std::list<Mesh::FaceHandle>& ungrouped_faces, TopologyGraph& graph);

	/// @brief Build a topology graph connecting all neighbouring regions.
	/// @param graph The graph to build.
	void buildTopologyGraph(TopologyGraph& graph);

	bool faceIsGrouped(const Mesh::FaceHandle fh, const TopologyGraph& graph) const;
	bool normalsAreCloseEnough(const Mesh::Normal& n_1, const Mesh::Normal& n_2, float threshold) const;

	/// Extends a given neighborhood with that of fh
	void extendNeighborhood(const Mesh::FaceHandle fh, std::list<Mesh::FaceHandle>& neighbors);
};

class EdgeCollapse {
	using RegionID = int;

private:
	OpenMesh::VPropHandleT<Quadric> v_quadric;
	Mesh* mesh{ nullptr };
	TopologyGraph* graph{ nullptr };
	std::map<RegionID, Quadric> region_quadrics;

	void computeRegionQuadrics();

	Quadric& vertexQuadric(const Mesh::VertexHandle vh) {
		return mesh->property(v_quadric, vh);
	}

public:
	EdgeCollapse(Mesh& mesh, TopologyGraph& graph) : mesh{ &mesh }, graph{ &graph } {
		this->mesh->add_property(v_quadric);
	}

	~EdgeCollapse() {
		this->mesh->remove_property(v_quadric);
	}

	Quadric computeQuadric(const Equation::Plane& plane);
	Quadric computeVertexQuadric(const Mesh::VertexHandle vh);
	Quadric computeEdgeQuadric(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1);

	void computeVerticesQuadrics();
	float computeEdgeError(const Vector3f& v, const Quadric& e_quadric);

	Vector3f computeNewVertex(const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1);
};