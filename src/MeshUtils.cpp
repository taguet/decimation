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


std::set<Mesh::VertexHandle> MeshUtils::getNeighboringVertices(const Mesh& mesh, const Mesh::VertexHandle vh) {
	std::set<Mesh::VertexHandle> neighbors{};
	for (auto& vv_it{ mesh.cvv_iter(vh) }; vv_it; ++vv_it) {
		neighbors.insert(vv_it);
	}
	return neighbors;
}


Mesh::HalfedgeHandle MeshUtils::findHalfedge(const Mesh& mesh, const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1)
{
	Mesh::HalfedgeHandle hh{ mesh.find_halfedge(vh_0, vh_1) };
	if (!hh.is_valid()) {
		hh = mesh.find_halfedge(vh_1, vh_0);
	}
	return hh;
}

Mesh::EdgeHandle MeshUtils::findEdge(const Mesh& mesh, const Mesh::VertexHandle vh_0, const Mesh::VertexHandle vh_1)
{
	Mesh::HalfedgeHandle hh{ findHalfedge(mesh, vh_0, vh_1) };

	if (!hh.is_valid())	
		return Mesh::InvalidEdgeHandle;

	return mesh.edge_handle(hh);
}


Mesh::FaceHandle MeshUtils::getAdjacentFace(const Mesh& mesh, const Mesh::FaceHandle fh, const Mesh::HalfedgeHandle hh) {
	assert(mesh.face_handle(hh) == fh);	// Make sure the halfedge belongs to the given face
	Mesh::HalfedgeHandle hh_opp{ mesh.opposite_halfedge_handle(hh) };
	return mesh.face_handle(hh_opp);
}


void MeshUtils::vertex_handles(const Mesh& mesh, const Mesh::EdgeHandle eh, Mesh::VertexHandle& vh_0, Mesh::VertexHandle& vh_1) {
	Mesh::HalfedgeHandle hh{ mesh.halfedge_handle(eh, 0) };
	if (!eh.is_valid())		hh = mesh.halfedge_handle(eh, 1);
	vh_0 = mesh.from_vertex_handle(hh);
	vh_1 = mesh.to_vertex_handle(hh);
}


void MeshUtils::projectVertexToLine(Mesh& mesh, const Mesh::VertexHandle vh, const Equation::Line& line) {
	const Mesh::Point p{ mesh.point(vh) };
	const Eigen::Vector3f projected{ line.projectPoint(Eigen::Vector3f{p[0], p[1], p[2]}) };
	mesh.set_point(vh, Mesh::Point{ projected[0], projected[1], projected[2] });
}


//Equation::Plane MeshUtils::fitPlaneToVertices(Mesh& mesh, std::set<Mesh::VertexHandle>& vertices) {
//	MatrixXf points{ vertices.size(), 3 }; // point coordinates
//	{
//		std::set<Mesh::VertexHandle>::iterator it{ vertices.begin() };
//		for (int i{ 0 }; i < points.rows(), it != vertices.end(); ++i, ++it) {
//			Mesh::Point p{ mesh.point(*it) };
//			points.row(i) = Vector3f{ p[0], p[1], p[2] };
//		}
//	}
//	MatrixXf XY{ 4,4 };
//	Vector4f Z{ 0, 0, 0, 0 };
//	XY.fill(0);
//
//	for (int i{ 0 }; i < points.rows(); ++i) {
//		float x{ points(i, 0) };
//		float y{ points(i, 1) };
//		float z{ points(i, 2) };
//		XY.row(0) += Vector4f{ 2*x*x, x*y, x*z, 2*x };
//		XY.row(1) += Vector4f{ x*y, 2*y*y, y*z, 2*y };
//		XY.row(2) += Vector4f{ x*z, y*z, 2*z*z, 2*z };
//		XY.row(3) += Vector4f{ 2*x, 2*y, 2*z, 3 };
//		Z += Vector4f{ x*y+x*z, x*y+y*z, x*z+y*z, x+y+z };
//	}
//	return Eigen::FullPivLU<MatrixXf>(XY).solve(Z);
//	//return  { Vector4f{XY.inverse() * Z} };
//}


Equation::Plane MeshUtils::fitPlaneToVertices(Mesh& mesh, std::set<Mesh::VertexHandle>& vertices) {
	// PCA
	MatrixXf points{ vertices.size(), 3 };
	{
		int p_i{ 0 };
		for (auto v_it{ vertices.begin() }; v_it != vertices.end(); ++v_it, ++p_i) {
			Mesh::Point p{ mesh.point(*v_it) };
			points.row(p_i) << p[0], p[1], p[2];
		}
	}

	MatrixXf centered{ points.rowwise() - points.colwise().mean() };	// Center the data by subtracting the mean of each column to each row

	////std::cerr << "means:" << points.colwise().mean() << "\n";
	MatrixXf cov{ (centered.transpose() * centered) / vertices.size() };
	Eigen::SelfAdjointEigenSolver<MatrixXf> eig(cov);
	MatrixXf eigenvectors{ eig.eigenvectors() };
	MatrixXf pca_transform{ eigenvectors.rightCols(3) };	// Last two columns are the direction vectors, third last is normal.
	// We can get the normal from the eigenvectors
	Vector3f normal{ pca_transform.col(0).normalized()};
	////std::cerr << "Normal:" << normal << "\n";
	// Get a point from the plane
	Vector3f point{ points.colwise().mean() };
	////std::cerr << "Point:" << point << "\n";
	return Equation::Plane(normal[0], normal[1], normal[2], -normal.dot(point));
}


Equation::Plane MeshUtils::computeFacePlane(const Mesh& mesh, const Mesh::FaceHandle fh) {
	Mesh::HalfedgeHandle heh{ mesh.halfedge_handle(fh) };
	Mesh::Normal vec_1, vec_2;
	mesh.calc_sector_vectors(heh, vec_1, vec_2);	//Compute AB and AC
	Mesh::Normal normal{ OpenMesh::cross(vec_1, vec_2) };
	normal.normalize();
	Mesh::Point p{ mesh.point(mesh.to_vertex_handle(heh)) }; //a, b and c

	float d{ normal[0] * p[0] + normal[1] * p[1] + normal[2] * p[2] };
	return Equation::Plane{ normal[0], normal[1], normal[2], -d };
}


Vector3f MeshUtils::toEigen(const Mesh::Point& p) {
	return Vector3f{ p[0], p[1], p[2] };
}