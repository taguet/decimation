#include "MeshTools.h"
#include <math.h>

void MeshTools::setMesh(Mesh& mesh) {
	this->mesh_ = mesh;
}

// Computes the area of a face from a given halfedge.
double MeshTools::computeFaceArea(const OpenMesh::HalfedgeHandle heh) {
	return mesh_.calc_sector_area(heh);
}

// Computes the vertex area of a given vertex, which corresponds to one third of the sum of all surrounding faces' areas.
double MeshTools::computeVertexArea(const OpenMesh::VertexHandle vh) {
	double area{ 0.0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		area += computeFaceArea(voh_it);
	}
	return area / 3.0;
}


double MeshTools::cotan(double angle) {
	return 1.0 / tan(angle);
}

// Computes the opposite angles to a given halfedge.
std::pair<double, double> MeshTools::getOppositeAngles(const OpenMesh::HalfedgeHandle oh) {
	OpenMesh::HalfedgeHandle opposite{ mesh_.opposite_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle left{ mesh_.next_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle right{ mesh_.next_halfedge_handle(opposite) };
	return { mesh_.calc_sector_angle(left), mesh_.calc_sector_angle(right) };
}

// Computes the discrete Laplacian of a vertex.
// function should be a scalar function.
template <class T>
T MeshTools::discreteLaplacian(const OpenMesh::VertexHandle vh, T(*function)(const OpenMesh::VertexHandle)) {
	T sum{};
	for (auto vv_it{ mesh_.vv_iter(vh) }; vv_it; ++vv_it) {
		std::pair<double, double> angles;
		sum += (cotan(angles.first) + cotan(angles.right)) * (function(vv_it) - function(vh));
	}
	return sum / (2 * computeVertexArea(vh));
}