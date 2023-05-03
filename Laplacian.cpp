#include "Laplacian.h"

Laplacian::Laplacian(Mesh& mesh) {
	this->mesh_ = &mesh;
}


Laplacian::~Laplacian() {
	delete this->mesh_;
}


//============== UNIFORM LAPLACIAN =================

OpenMesh::Vec3f UniformLaplacian::computeLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	Mesh::Point p_i{ mesh_->point(vh) };
	int i{ 0 };
	for (auto voh_it{ mesh_->voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		Mesh::Point p_neighbour{ mesh_->point(mesh_->to_vertex_handle(voh_it)) };
		sum += p_neighbour - p_i;
	}
	if (i == 0)
		return sum;
	else
		return sum / i;
}


//============== COTANGENT LAPLACIAN =================

OpenMesh::Vec3f CotangentLaplacian::computeLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	Mesh::Point p_i{ mesh_->point(vh) };
	for (auto voh_it{ mesh_->voh_iter(vh) }; voh_it; ++voh_it) {
		Mesh::Point p_neighbour{ mesh_->point(mesh_->to_vertex_handle(voh_it)) };
		std::pair<float, float> angles{ getOppositeAngles(voh_it) };
		float weight{ cotan(angles.first) + cotan(angles.second) };
		OpenMesh::Vec3f vec{ p_neighbour - p_i };
		sum += weight * vec;
	}
	return sum / (2.0f * computeVertexArea(vh));
}


//============== ANISOTROPIC LAPLACIAN =================

OpenMesh::Vec3f AnisotropicLaplacian::computeLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	Mesh::Point p_i{ mesh_->point(vh) };
	int i{ 0 };
	for (auto voh_it{ mesh_->voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		Mesh::Point centroid{ computeCentroid(voh_it) };
		OpenMesh::FaceHandle face{ mesh_->face_handle(voh_it) };
		if (!face.is_valid())	continue;
		OpenMesh::Vec3f f_normal{ filterFaceNormal(face) };
		OpenMesh::Vec3f normal{ mesh_->normal(face) };
		OpenMesh::Vec3f toCentroid{ centroid - p_i };
		sum += dot(toCentroid, f_normal) * normal;
	}
	if (i == 0)
		return sum;
	else
		return sum / i;
}
