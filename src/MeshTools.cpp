#include "MeshTools.h"
#include <math.h>


MeshTools::MeshTools(void) {
	
}

void MeshTools::setMesh(Mesh& mesh) {
	this->mesh_ = mesh;
	this->mesh_.add_property(laplacian);
}

const Mesh* MeshTools::getMesh() {
	return &mesh_;
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
	double area{ 0.0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		area += computeFaceArea(voh_it);
	}
	return area / 3.0;
}


float MeshTools::cotan(double angle) {
	return 1.0 / tan(angle);
}

// Computes the opposite angles to a given halfedge.
std::pair<float, float> MeshTools::getOppositeAngles(const OpenMesh::HalfedgeHandle oh) {
	OpenMesh::HalfedgeHandle opposite{ mesh_.opposite_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle left{ mesh_.next_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle right{ mesh_.next_halfedge_handle(opposite) };
	return { mesh_.calc_sector_angle(left), mesh_.calc_sector_angle(right) };
}

OpenMesh::Vec3f MeshTools::discreteLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{0.0f, 0.0f, 0.0f};
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		std::pair<float, float> angles{getOppositeAngles(voh_it)};
		sum += (cotan(angles.first) + cotan(angles.second)) * (mesh_.point(mesh_.to_vertex_handle(voh_it)) - mesh_.point(vh));
	}
	return sum / (2.0 * computeVertexArea(vh));
}


OpenMesh::Vec3f MeshTools::uniformLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	int i{ 0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		sum += mesh_.point(mesh_.to_vertex_handle(voh_it)) - mesh_.point(vh);
	}
	return sum / i;
}


void MeshTools::smoothMesh(int h, float lambda) {
	for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
		auto p{ mesh_.point(v_it) };
		mesh_.set_point(v_it, position(v_it) + laplacian_displacement(v_it) * h * lambda);
	}
}

void MeshTools::smoothMesh(int iterations=1) {
	for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
		OpenMesh::Vec3f lapl{ laplacian_displacement(v_it) };
		for (int i{ 1 }; i < iterations; ++i) {
			lapl *= laplacian_displacement(v_it);
		}
		auto p{ mesh_.point(v_it) };
		mesh_.set_point(v_it, p + lapl);
	}
}