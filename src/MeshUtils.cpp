#include "MeshUtils.h"

// Computes the area of a face from a given halfedge.
float MeshUtils::computeFaceArea(Mesh& mesh, const OpenMesh::HalfedgeHandle heh) {
	return mesh.calc_sector_area(heh);
}

float MeshUtils::computeFaceArea(Mesh& mesh, const OpenMesh::FaceHandle fh) {
	Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(fh) };
	return computeFaceArea(mesh, heh);
}

// Computes the vertex area of a given vertex, which corresponds to one third of the sum of all surrounding faces' areas.
float MeshUtils::computeVertexArea(Mesh& mesh, const OpenMesh::VertexHandle vh) {
	float area{ 0.0f };
	for (auto voh_it{ mesh.voh_iter(vh) }; voh_it; ++voh_it) {
		area += computeFaceArea(mesh, voh_it);
	}
	return area / 3.0f;
}


// Computes the opposite angles to a given halfedge.
std::pair<float, float> MeshUtils::getOppositeAngles(Mesh& mesh, const OpenMesh::HalfedgeHandle oh) {
	OpenMesh::HalfedgeHandle opposite{ mesh.opposite_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle left{ mesh.next_halfedge_handle(oh) };
	OpenMesh::HalfedgeHandle right{ mesh.next_halfedge_handle(opposite) };
	return { computeAngle(mesh, left), computeAngle(mesh, right) };
}


float MeshUtils::computeAngle(Mesh& mesh, const OpenMesh::HalfedgeHandle heh) {
	OpenMesh::Vec3f e1, e2;
	mesh.calc_sector_vectors(heh, e1, e2);
	return acosf(dot(e1, e2) / (e1.norm() * e2.norm()));
}


float MeshUtils::cotan(float angle) {
	return 1.0f / tanf(angle);
}


OpenMesh::Vec3f MeshUtils::computeCentroid(Mesh& mesh, const OpenMesh::HalfedgeHandle heh) {
	OpenMesh::VertexHandle v0, v1, v2;
	v0 = mesh.from_vertex_handle(heh);
	v1 = mesh.to_vertex_handle(heh);
	v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
	Mesh::Point p0, p1, p2;
	p0 = mesh.point(v0);
	p1 = mesh.point(v1);
	p2 = mesh.point(v2);
	return (p0 + p1 + p2) / 3.0f;
}


float MeshUtils::computeAverageEdgeLength(Mesh& mesh) {
	float avrg_edge_length{ 0.0f };
	for (auto e_it{ mesh.edges_begin() }; e_it != mesh.edges_end(); ++e_it) {
		avrg_edge_length += mesh.calc_edge_length(e_it);
	}
	return avrg_edge_length / mesh.n_edges();
}


float MeshUtils::computeVarianceEdgeLength(Mesh& mesh) {
	float avrg_edge_length{ computeAverageEdgeLength(mesh) };
	float variance{ 0.0f };
	for (auto e_it{ mesh.edges_begin() }; e_it != mesh.edges_end(); ++e_it) {
		variance += powf(mesh.calc_edge_length(e_it) - avrg_edge_length, 2);
	}
	return variance / mesh.n_edges();
}


void MeshUtils::getFaceNeighbors(Mesh& mesh, Mesh::FaceHandle fh, std::list<Mesh::FaceHandle>& neighbors) {
	for (auto f_iter{ mesh.ff_iter(fh) }; f_iter; ++f_iter) {
		neighbors.push_back(f_iter);
	}
}


Vector3f MeshUtils::fitPlaneToVertices(Mesh& mesh, std::set<Mesh::VertexHandle>& vertices) {
	MatrixXf regressors{ vertices.size(), 3 };	// x & y
	VectorXf observed{ vertices.size() };	// z
	Vector3f parameters{ 0.0f, 0.0f, 0.0f };	//c, a & b

	{
		std::set<Mesh::VertexHandle>::iterator it{ vertices.begin() };
		for (int i{ 0 }; i < regressors.rows(), it != vertices.end(); ++i, ++it) {
			Mesh::Point p{ mesh.point(*it) };
			regressors(i, 0) = 1.0f;
			regressors(i, 1) = p[0];
			regressors(i, 2) = p[1];
			observed(i) = p[2];
		}
	}
	MatrixXf t_regressors{ regressors.transpose() };
	parameters = (t_regressors * regressors).inverse() * t_regressors * observed;
	return parameters;
}
