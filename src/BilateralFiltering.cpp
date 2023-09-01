#include "BilateralFiltering.h"

BilateralFiltering::BilateralFiltering(Mesh& mesh) : mesh{ &mesh } {
	variance_edge_length = MeshUtils::computeVarianceEdgeLength(mesh);
	this->mesh->add_property(filtered_normals);
}


BilateralFiltering::~BilateralFiltering() {
	mesh->remove_property(filtered_normals);
}


void BilateralFiltering::initializeNormals() {
	for (auto fh_it{ mesh->faces_begin() }; fh_it != mesh->faces_end(); ++fh_it) {
		filteredNormal(fh_it) = mesh->normal(fh_it) / mesh->normal(fh_it).norm();
	}
}


void BilateralFiltering::filterFaceNormals(int iterations, float threshold) {
	initializeNormals();
	std::vector<Mesh::Normal> f_normals;	// Temporary container to store filtered normals as they are computed
	f_normals.reserve(mesh->n_faces());
	for (int i{ 0 }; i < iterations; ++i) {
		for (auto fh_it{ mesh->faces_begin() }; fh_it != mesh->faces_end(); ++fh_it) {
			f_normals.push_back(filterFaceNormal(fh_it, threshold));
		}
		for (int f_id{ 0 }; f_id < mesh->n_faces(); ++f_id) {	// Apply filtered normals
			filteredNormal(mesh->face_handle(f_id)) = f_normals.at(f_id);
		}
	}
}


Mesh::Normal BilateralFiltering::filterFaceNormal(const OpenMesh::FaceHandle fh, float threshold) {
	float area{ MeshUtils::computeFaceArea(*mesh, mesh->halfedge_handle(fh)) };
	Mesh::Normal f_normal{ weighFaceNormal(fh, fh, area, threshold) };
	for (auto ff_it{ mesh->ff_iter(fh) }; ff_it; ++ff_it) {
		if (dot(filteredNormal(ff_it), filteredNormal(fh)) < cos(threshold))
			continue;
		f_normal += weighFaceNormal(fh, ff_it, area, threshold);
	}
	return f_normal / f_normal.norm();
}


float BilateralFiltering::computeDistanceWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2) {
	Mesh::Point c_0, c_1;
	c_0 = MeshUtils::computeCentroid(*mesh, mesh->halfedge_handle(fh1));
	c_1 = MeshUtils::computeCentroid(*mesh, mesh->halfedge_handle(fh2));
	float sqr_dst{ (c_0 - c_1).sqrnorm() };
	return exp(-sqr_dst / (2 * variance_edge_length));
}


float BilateralFiltering::computeProximityWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2, float threshold) {
	Mesh::Normal n_0, n_1;
	n_0 = filteredNormal(fh1);
	n_1 = filteredNormal(fh2);
	float num{ pow(1 - dot(n_0, n_1), 2.0f) };
	float denom{ pow(1 - cos(threshold), 2.0f) };
	return exp(-num / denom);
}


void BilateralFiltering::update() {
	variance_edge_length = MeshUtils::computeVarianceEdgeLength(*mesh);
}


Mesh::Normal& BilateralFiltering::filteredNormal(const Mesh::FaceHandle fh) {
	return mesh->property(filtered_normals, fh);
}


Mesh::Normal BilateralFiltering::weighFaceNormal(const OpenMesh::FaceHandle fh_i, const OpenMesh::FaceHandle fh_j, const float area, const float threshold) {
	float alpha{ computeDistanceWeight(fh_i, fh_j) };
	float beta{ computeProximityWeight(fh_i, fh_j, threshold) };
	return area * alpha * beta * filteredNormal(fh_j);
}