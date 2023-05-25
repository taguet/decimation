#pragma once

#include "gl.h"
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <map>
#include <set>
#include "Laplacian.h"
#include <Eigen/Dense>


typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::Vector3f;


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
			std::cout << "Iteration " << i + 1 << std::endl;
			laplacian->computeLaplacians();
			for (auto v_it{ mesh_->vertices_begin() }; v_it != mesh_->vertices_end(); ++v_it) {
				mesh_->set_point(v_it, mesh_->point(v_it) + laplacian->laplacian_displacement(v_it) * factor);
			}
			mesh_->update_normals();
		}
	}


	void extractRegions();


	int& faceGroup(Mesh::FaceHandle fh) {
		return mesh_->property(f_group, fh);
	}


	Vector3f fitPlaneToVertices(std::set<Mesh::VertexHandle>& vertices) const;


private:
	class TopologyGraph {
	private:
		class Node {
		public:
			const int id;

			Node(TopologyGraph& parent, int id) : id{ id }, parent{ &parent } {}
			void add(Mesh::FaceHandle fh);
			void regroupIntoSelf(Node& region);
			bool empty() const { return faces.empty(); }
			float computeArea() const;
			void fitPlane();
			float sumVertexProjectedDistances();
			void updateIndices();	//TODO delete
		private:
			TopologyGraph* parent{ nullptr };
			std::set<Mesh::FaceHandle> faces;
			std::set<Mesh::VertexHandle> vertices;
			Vector3f plane_params;
		};

		std::map<int, Node> regions;
		std::map<int, std::set<int>> edges;
		MeshTools* parent{ nullptr };

		int findTargetRegion(int regionID, float fitting_threshold);
		void regroupRegionIntoTarget(int regionID, int targetID);
		void ungroupRegion(int regionID);

	public:
		const float area_threshold;
		const float fitting_threshold;

		TopologyGraph(MeshTools& parent, float area_threshold, float fitting_threshold);
		int size() {
			return regions.size();
		}

		void addFaceToRegion(int regionID, Mesh::FaceHandle fh);

		Node& getRegion(int regionID) { return regions.at(regionID); }

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
		void fitPlanes();
		void updateIndices();	//TODO delete
	};

	OpenMesh::FPropHandleT<int> f_group;

	/// @brief Region growing algorithm to detect and separate all planar regions.
	/// @param ungrouped_faces Every faces in the mesh that haven't been sorted in a group yet.
	/// @param graph The topology graph to update throughout the process.
	void growRegions(std::list<Mesh::FaceHandle>& ungrouped_faces, TopologyGraph& graph);

	/// @brief Build a topology graph connecting all neighbouring regions.
	/// @param graph The graph to build.
	void buildTopologyGraph(TopologyGraph& graph);
};

