#include "MeshTools.h"
#include <math.h>


MeshTools::MeshTools(void) {
	
}

void MeshTools::setMesh(Mesh& mesh) {
	this->mesh_ = mesh;
}

OpenMesh::Vec3f MeshTools::position(const OpenMesh::VertexHandle vh) {
	auto p{ mesh_.point(vh) };
	return {p[0], p[1], p[2]};
}

// Computes the area of a face from a given halfedge.
double MeshTools::computeFaceArea(const OpenMesh::HalfedgeHandle heh) {
	return mesh_.calc_sector_area(heh);
}

// Computes the vertex area of a given vertex, which corresponds to one third of the sum of all surrounding faces' areas.
double MeshTools::computeVertexArea(const OpenMesh::VertexHandle vh) {
	std::cout << position(vh);
	std::cout << vh.is_valid() << std::endl;
	double area{ 0.0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		area += computeFaceArea(voh_it);
	}
	return area / 3.0;
}


double MeshTools::cotan(double angle) {
	return 1.0 / tan(angle);
}

//// Computes the opposite angles to a given halfedge.
std::pair<float, float> MeshTools::getOppositeAngles(const OpenMesh::HalfedgeHandle oh) {
	OpenMesh::HalfedgeHandle opposite{ mesh_.opposite_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle left{ mesh_.next_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle right{ mesh_.next_halfedge_handle(opposite) };
	return { mesh_.calc_sector_angle(left), mesh_.calc_sector_angle(right) };
}

OpenMesh::Vec3f MeshTools::discreteLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{0.0f, 0.0f, 0.0f};
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		//std::pair<float, float> angles{getOppositeAngles(voh_it)};
		sum += /*(cotan(angles.first) + cotan(angles.right)) * */(position(mesh_.to_vertex_handle(voh_it)) - position(vh));
	}
	return sum / (2 * computeVertexArea(vh));
}