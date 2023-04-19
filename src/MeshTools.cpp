#include "MeshTools.h"
#include <math.h>


MeshTools::MeshTools(void) {
	
}

void MeshTools::setMesh(Mesh& mesh) {
	this->mesh_ = mesh;
	this->mesh_.add_property(laplacian);
	this->variance_edge_length = computeVarianceEdgeLength();
}


float MeshTools::computeAverageEdgeLength() {
	float avrg_edge_length{ 0.0f };
	for (auto e_it{ mesh_.edges_begin() }; e_it != mesh_.edges_end(); ++e_it) {
		avrg_edge_length += mesh_.calc_edge_length(e_it);
	}
	return avrg_edge_length / mesh_.n_edges();
}


float MeshTools::computeVarianceEdgeLength() {
	float avrg_edge_length{ computeAverageEdgeLength() };
	float variance{ 0.0f };
	for (auto e_it{ mesh_.edges_begin() }; e_it != mesh_.edges_end(); ++e_it) {
		variance += pow(mesh_.calc_edge_length(e_it) - avrg_edge_length, 2);
	}
	return variance / mesh_.n_edges();
}


const Mesh* MeshTools::getMesh() {
	return &mesh_;
}

OpenMesh::Vec3f MeshTools::position(const OpenMesh::VertexHandle vh) {
	auto p{ mesh_.point(vh) };
	return {p[0], p[1], p[2]};
}

// Computes the area of a face from a given halfedge.
float MeshTools::computeFaceArea(const OpenMesh::HalfedgeHandle heh) {
	return mesh_.calc_sector_area(heh);
}

// Computes the vertex area of a given vertex, which corresponds to one third of the sum of all surrounding faces' areas.
float MeshTools::computeVertexArea(const OpenMesh::VertexHandle vh) {
	double area{ 0.0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		area += computeFaceArea(voh_it);
	}
	return area / 3.0;
}


float MeshTools::cotan(double angle) {
	if (angle == 0)	return 0;
	return 1.0 / tan(angle);
}

// Computes the opposite angles to a given halfedge.
std::pair<float, float> MeshTools::getOppositeAngles(const OpenMesh::HalfedgeHandle oh) {
	OpenMesh::HalfedgeHandle opposite{ mesh_.opposite_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle left{ mesh_.next_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle right{ mesh_.next_halfedge_handle(opposite) };
	return { computeAngle(left), computeAngle(right) };
}


float MeshTools::computeAngle(const OpenMesh::HalfedgeHandle heh) {
	OpenMesh::Vec3f e1, e2;
	mesh_.calc_sector_vectors(heh, e1, e2);
	return cos(OpenMesh::dot(e1.normalize(), e2.normalize()));
}


OpenMesh::Vec3f MeshTools::computeCentroid(const OpenMesh::HalfedgeHandle heh) {
	OpenMesh::VertexHandle v0, v1, v2;
	v0 = mesh_.from_vertex_handle(heh);
	v1 = mesh_.to_vertex_handle(heh);
	v2 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heh));
	Mesh::Point p0, p1, p2;
	p0 = mesh_.point(v0);
	p1 = mesh_.point(v1);
	p2 = mesh_.point(v2);
	return (p0 + p1 + p2) / 3.0f;
}


Mesh::Normal MeshTools::filterFaceNormal(const OpenMesh::FaceHandle fh) {
	float area{ computeFaceArea(mesh_.halfedge_handle(fh)) };
	Mesh::Normal normal{ mesh_.normal(fh) };
	//TODO weights
}


float MeshTools::computeDistanceWeight(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle neighbour) {
	Mesh::Point c_0, c_1;
	c_0 = computeCentroid(mesh_.halfedge_handle(fh));
	c_1 = computeCentroid(mesh_.halfedge_handle(neighbour));
	float sqr_dst{ (c_0 - c_1).sqrnorm() };
	return exp(-sqr_dst / (2 * variance_edge_length));
}


OpenMesh::Vec3f MeshTools::cotangentLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{0.0f, 0.0f, 0.0f};
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		std::pair<float, float> angles{getOppositeAngles(voh_it)};
		float weight{ (cotan(angles.first) + cotan(angles.second)) };
		OpenMesh::Vec3f vec{ (mesh_.point(mesh_.to_vertex_handle(voh_it)) - mesh_.point(vh)) };
		sum += weight * vec;
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


OpenMesh::Vec3f MeshTools::anisotropicLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	int i{ 0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		Mesh::Point centroid{ computeCentroid(voh_it) };
		//TODO: Filtered normal
		Mesh::Normal n{ mesh_.normal(mesh_.face_handle(voh_it)) };
	}
	return sum; //TODO
}


void MeshTools::taubinSmoothing(float lambda, float mu) {
	for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
		mesh_.set_point(v_it, position(v_it) + laplacian_displacement(v_it) * lambda);
		mesh_.set_point(v_it, position(v_it) + laplacian_displacement(v_it) * mu);
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