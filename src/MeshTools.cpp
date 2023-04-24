#include "MeshTools.h"
#include <math.h>


MeshTools::MeshTools(void) {
	
}

void MeshTools::setMesh(Mesh& mesh) {
	this->mesh_ = mesh;
	this->mesh_.add_property(laplacian);
}


void MeshTools::calc_discrete_laplacian() {
	for (auto v_iter{ mesh_.vertices_begin() }; v_iter != mesh_.vertices_end(); ++v_iter) {
		//laplacian_displacement(v_iter) = anisotropicLaplacian(v_iter);
		laplacian_displacement(v_iter) = cotangentLaplacian(v_iter);
		//laplacian_displacement(v_iter) = uniformLaplacian(v_iter);
		//std::cout << laplacian_displacement(v_iter) << std::endl;
	}
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
		variance += powf(mesh_.calc_edge_length(e_it) - avrg_edge_length, 2);
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
	float area{ 0.0f };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		area += computeFaceArea(voh_it);
	}
	return area / 3.0f;
}


float MeshTools::cotan(float angle) {
	return 1.0f / tanf(angle);
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
	return acosf(dot(e1.normalize(), e2.normalize()));
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


Mesh::Normal MeshTools::filterFaceNormal(const OpenMesh::FaceHandle fh, float threshold) {
	float area{ computeFaceArea(mesh_.halfedge_handle(fh)) };
	Mesh::Normal f_normal{weighFaceNormal(fh, fh, area)};
	for (auto ff_it{ mesh_.ff_iter(fh) }; ff_it; ++ff_it) {
		if (dot(mesh_.normal(ff_it), mesh_.normal(fh)) >= cos(threshold))
			continue;
		f_normal += weighFaceNormal(fh, ff_it, area);
	}
	return f_normal / f_normal.norm();
}


Mesh::Normal MeshTools::weighFaceNormal(const OpenMesh::FaceHandle fh_i, const OpenMesh::FaceHandle fh_j, const float area) {
	float alpha{ computeDistanceWeight(fh_i, fh_j) };
	float beta{ computeProximityWeight(fh_i, fh_j) };
	return area * alpha * beta * mesh_.normal(fh_j);
}


float MeshTools::computeDistanceWeight(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle neighbour) {
	Mesh::Point c_0, c_1;
	c_0 = computeCentroid(mesh_.halfedge_handle(fh));
	c_1 = computeCentroid(mesh_.halfedge_handle(neighbour));
	float sqr_dst{ (c_0 - c_1).sqrnorm() };
	return exp(-sqr_dst / (2 * variance_edge_length));
}


float MeshTools::computeProximityWeight(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle neighbour, float threshold) {
	Mesh::Normal n_0, n_1;
	n_0 = mesh_.normal(fh);
	n_1 = mesh_.normal(neighbour);
	float num{pow(1 - dot(n_0, n_1), 2.0f)};
	float denom{pow(1 - cos(threshold), 2.0f)};
	return exp( - num / denom );
}


OpenMesh::Vec3f MeshTools::cotangentLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{0.0f, 0.0f, 0.0f};
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it) {
		std::pair<float, float> angles{getOppositeAngles(voh_it)};
		float weight{ cotan(angles.first) + cotan(angles.second) };
		OpenMesh::Vec3f vec{ mesh_.point(mesh_.to_vertex_handle(voh_it)) - mesh_.point(vh) };
		//std::cout << "sum=" << sum << "\tvec=" << vec << "\tweight=" << weight << std::endl;
		sum += weight * vec;
	}
	//std::cout << "sum=" << sum << "\tarea=" << computeVertexArea(vh) << "\tL=" << sum / (2.0f * computeVertexArea(vh)) << std::endl;
	return sum / (2.0f * computeVertexArea(vh));
}


OpenMesh::Vec3f MeshTools::uniformLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	int i{ 0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		sum += mesh_.point(mesh_.to_vertex_handle(voh_it)) - mesh_.point(vh);
	}
	if (i == 0)
		return sum;
	else
		return sum / i;
}


OpenMesh::Vec3f MeshTools::anisotropicLaplacian(const OpenMesh::VertexHandle vh) {
	OpenMesh::Vec3f sum{ 0.0f, 0.0f, 0.0f };
	int i{ 0 };
	for (auto voh_it{ mesh_.voh_iter(vh) }; voh_it; ++voh_it, ++i) {
		//std::cout << "Vi=" << vh.idx() << "\tVj=" << mesh_.to_vertex_handle(voh_it).idx() << std::endl;
		Mesh::Point centroid{ computeCentroid(voh_it) };
		OpenMesh::Vec3f f_normal{ filterFaceNormal(mesh_.face_handle(voh_it)) };
		OpenMesh::Vec3f normal{ mesh_.normal(mesh_.face_handle(voh_it)) };
		OpenMesh::Vec3f toCentroid{ centroid - mesh_.point(vh) };
		sum += dot(toCentroid, f_normal) * normal;
	}
	if (i == 0)	
		return sum;
	else
		return sum / i;
}


void MeshTools::taubinSmoothing(float lambda, float mu, int iterations) {
	for (int i{ 0 }; i < iterations; ++i) {
		variance_edge_length = computeVarianceEdgeLength();
		calc_discrete_laplacian();
		for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
			mesh_.set_point(v_it, mesh_.point(v_it) + laplacian_displacement(v_it) * lambda);
		}
		variance_edge_length = computeVarianceEdgeLength();
		mesh_.update_normals();
		calc_discrete_laplacian();
		for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
			mesh_.set_point(v_it, mesh_.point(v_it) + laplacian_displacement(v_it) * mu);
		}
		mesh_.update_normals();
	}
}

void MeshTools::smoothMesh(int iterations) {
	for (int i{ 0 }; i < iterations; ++i) {
		std::cout << "Iteration " << i + 1 << std::endl;
		variance_edge_length = computeVarianceEdgeLength();
		calc_discrete_laplacian();
		for (auto v_it{ mesh_.vertices_begin() }; v_it != mesh_.vertices_end(); ++v_it) {
			//std::cout << "Point: " << mesh_.point(v_it) << "\tLaplacien: " << laplacian_displacement(v_it) << std::endl;
			mesh_.set_point(v_it, mesh_.point(v_it) + laplacian_displacement(v_it));
		}
		mesh_.update_normals();
	}
}