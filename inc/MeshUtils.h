#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class MeshUtils
{
public:
	static float computeFaceArea(Mesh& mesh, const OpenMesh::HalfedgeHandle heh);
	static float computeVertexArea(Mesh& mesh, const OpenMesh::VertexHandle vh);

	static std::pair<float, float> getOppositeAngles(Mesh& mesh, const OpenMesh::HalfedgeHandle oh);
	static float computeAngle(Mesh& mesh, const OpenMesh::HalfedgeHandle heh);
	static float cotan(float angle);

	static OpenMesh::Vec3f computeCentroid(Mesh& mesh, const OpenMesh::HalfedgeHandle heh);

	static float computeAverageEdgeLength(Mesh& mesh);
	static float computeVarianceEdgeLength(Mesh& mesh);
};

