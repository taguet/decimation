#include "MeshUtils.h"

// Computes the area of a face from a given halfedge.
float MeshUtils::computeFaceArea(Mesh& mesh, const OpenMesh::HalfedgeHandle heh) {
	return mesh.calc_sector_area(heh);
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