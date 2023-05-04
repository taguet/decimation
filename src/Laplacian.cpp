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
	return sum / (2.0f * MeshUtils::computeVertexArea(*mesh_, vh));
}


//============== ANISOTROPIC LAPLACIAN =================

AnisotropicLaplacian::AnisotropicLaplacian(Mesh& mesh) : Laplacian(mesh) {
	variance_edge_length = MeshUtils::computeVarianceEdgeLength(*mesh_);
}


OpenMesh::Vec3f AnisotropicLaplacian::computeLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	Mesh::Point p_i{ mesh_->point(vh) };
	int i{ 0 };
	for (auto voh_it{ mesh_->voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		Mesh::Point centroid{ MeshUtils::computeCentroid(*mesh_, voh_it) };
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


void AnisotropicLaplacian::computeLaplacians() {
	Laplacian::computeLaplacians();
	variance_edge_length = MeshUtils::computeVarianceEdgeLength(*mesh_);
}


Mesh::Normal AnisotropicLaplacian::filterFaceNormal(const OpenMesh::FaceHandle fh, float threshold) {
	float area{ MeshUtils::computeFaceArea(*mesh_, mesh_->halfedge_handle(fh)) };
	Mesh::Normal f_normal{ weighFaceNormal(fh, fh, area) };
	for (auto ff_it{ mesh_->ff_iter(fh) }; ff_it; ++ff_it) {
		if (dot(mesh_->normal(ff_it), mesh_->normal(fh)) < cos(threshold))
			continue;
		f_normal += weighFaceNormal(fh, ff_it, area);
	}
	return f_normal / f_normal.norm();
}


float AnisotropicLaplacian::computeDistanceWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2) {
	Mesh::Point c_0, c_1;
	c_0 = MeshUtils::computeCentroid(*mesh_, mesh_->halfedge_handle(fh1));
	c_1 = MeshUtils::computeCentroid(*mesh_, mesh_->halfedge_handle(fh2));
	float sqr_dst{ (c_0 - c_1).sqrnorm() };
	return exp(-sqr_dst / (2 * variance_edge_length));
}


float AnisotropicLaplacian::computeProximityWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2, float threshold) {
	Mesh::Normal n_0, n_1;
	n_0 = mesh_->normal(fh1);
	n_1 = mesh_->normal(fh2);
	float num{ pow(1 - dot(n_0, n_1), 2.0f) };
	float denom{ pow(1 - cos(threshold), 2.0f) };
	return exp(-num / denom);
}


Mesh::Normal AnisotropicLaplacian::weighFaceNormal(const OpenMesh::FaceHandle fh_i, const OpenMesh::FaceHandle fh_j, const float area) {
	float alpha{ computeDistanceWeight(fh_i, fh_j) };
	float beta{ computeProximityWeight(fh_i, fh_j) };
	return area * alpha * beta * mesh_->normal(fh_j);
}