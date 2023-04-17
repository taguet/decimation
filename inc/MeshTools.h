#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>


typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class MeshTools
{
private:
	Mesh mesh_;
public:
	
	MeshTools(void);
	void setMesh(Mesh &mesh);
	OpenMesh::Vec3f position(const OpenMesh::VertexHandle vh);
	double computeFaceArea(const OpenMesh::HalfedgeHandle heh);
	double computeVertexArea(const OpenMesh::VertexHandle vh);
	double cotan(double angle);
	std::pair<float, float> getOppositeAngles(const OpenMesh::HalfedgeHandle oh);

	OpenMesh::Vec3f discreteLaplacian(const OpenMesh::VertexHandle vh);
	
	// Computes the discrete Laplacian of a vertex.
	// function should be a scalar function.
	//template <typename T, typename F>
	//T discreteLaplacian(const OpenMesh::VertexHandle vh, T (F::*function)(const OpenMesh::VertexHandle)) {
	//	OpenMesh::Vec3f sum{};
	//	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
	//		//std::pair<double, double> angles{getOppositeAngles(voh_it)};
	//		sum += /*(cotan(angles.first) + cotan(angles.right)) * */(function(mesh_.to_vertex_handle(voh_it)) - function(vh));
	//	}
	//	return sum / (2 * computeVertexArea(vh));
	//}
};

