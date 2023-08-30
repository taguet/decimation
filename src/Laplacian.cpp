#include "Laplacian.h"

Laplacian::Laplacian(Mesh& mesh) {
	this->mesh_ = &mesh;
	this->mesh_->add_property(laplacian);
}


OpenMesh::Vec3f& Laplacian::laplacian_displacement(Mesh::VertexHandle _vh)
{
	return mesh_->property(laplacian, _vh);
}


void Laplacian::computeLaplacians() {
	for (auto v_iter{ mesh_->vertices_begin() }; v_iter != mesh_->vertices_end(); ++v_iter) {
		laplacian_displacement(v_iter) = computeLaplacian(v_iter);
	}
}


//============== UNIFORM LAPLACIAN =================

UniformLaplacian::UniformLaplacian(Mesh& mesh) : Laplacian(mesh) {}


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

CotangentLaplacian::CotangentLaplacian(Mesh& mesh) : Laplacian(mesh) {}


OpenMesh::Vec3f CotangentLaplacian::computeLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	Mesh::Point p_i{ mesh_->point(vh) };
	for (auto voh_it{ mesh_->voh_iter(vh) }; voh_it; ++voh_it) {
		Mesh::Point p_neighbour{ mesh_->point(mesh_->to_vertex_handle(voh_it)) };
		std::pair<float, float> angles{ MeshUtils::getOppositeAngles(*mesh_, voh_it) };
		float weight{ MeshUtils::cotan(angles.first) + MeshUtils::cotan(angles.second) };
		OpenMesh::Vec3f vec{ p_neighbour - p_i };
		sum += weight * vec;
	}
	return (sum / (2.0f * MeshUtils::computeVertexArea(*mesh_, vh))).normalize();
}


//============== ANISOTROPIC LAPLACIAN =================

AnisotropicLaplacian::AnisotropicLaplacian(Mesh& mesh) : filter{mesh}, Laplacian(mesh) {
}


OpenMesh::Vec3f AnisotropicLaplacian::computeLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	Mesh::Point p_i{ mesh_->point(vh) };
	int i{ 0 };
	for (auto voh_it{ mesh_->voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		Mesh::Point centroid{ MeshUtils::computeCentroid(*mesh_, voh_it) };
		OpenMesh::FaceHandle face{ mesh_->face_handle(voh_it) };
		if (!face.is_valid()) {
			--i;
			continue;
		}
		OpenMesh::Vec3f f_normal{ filter.filterFaceNormal(face) };
		OpenMesh::Vec3f normal{ mesh_->normal(face) };
		OpenMesh::Vec3f toCentroid{ centroid - p_i };
		sum += dot(toCentroid, f_normal) * normal;
	}
	if (i <= 0)
		return sum;
	else
		return sum / i;
}


void AnisotropicLaplacian::computeLaplacians() {
	Laplacian::computeLaplacians();
	filter.update();
}
