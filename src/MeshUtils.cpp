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


void MeshUtils::projectVertexToLine(Mesh& mesh, const Mesh::VertexHandle vh, const Equation::Line& line) {
	const Mesh::Point p{ mesh.point(vh) };
	const Eigen::Vector3f projected{ line.projectPoint(Eigen::Vector3f{p[0], p[1], p[2]}) };
	mesh.set_point(vh, Mesh::Point{ projected[0], projected[1], projected[2] });
}


Equation::Plane MeshUtils::fitPlaneToVertices(Mesh& mesh, std::set<Mesh::VertexHandle>& vertices) {
	MatrixXf points{ vertices.size(), 3 }; // point coordinates
	{
		std::set<Mesh::VertexHandle>::iterator it{ vertices.begin() };
		for (int i{ 0 }; i < points.rows(), it != vertices.end(); ++i, ++it) {
			Mesh::Point p{ mesh.point(*it) };
			points.row(i) = Vector3f{ p[0], p[1], p[2] };
		}
	}
	MatrixXf XY{ 4,4 };
	Vector4f Z{ 0, 0, 0, 0 };
	XY.fill(0);

	for (int i{ 0 }; i < points.rows(); ++i) {
		float x{ points(i, 0) };
		float y{ points(i, 1) };
		float z{ points(i, 2) };
		XY.row(0) += Vector4f{ 2*x*x, x*y, x*z, 2*x };
		XY.row(1) += Vector4f{ x*y, 2*y*y, y*z, 2*y };
		XY.row(2) += Vector4f{ x*z, y*z, 2*z*z, 2*z };
		XY.row(3) += Vector4f{ 2*x, 2*y, 2*z, 3 };
		Z += Vector4f{ x*y+x*z, x*y+y*z, x*z+y*z, x+y+z };
	}
	return  { Vector4f{XY.inverse() * Z} };
}
