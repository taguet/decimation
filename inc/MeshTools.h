#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>


typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class MeshTools
{
private:
	Mesh mesh_;
public:
	
	MeshTools(void);
	void setMesh(Mesh &mesh);
	const Mesh* getMesh();
	OpenMesh::Vec3f position(const OpenMesh::VertexHandle vh);
	float computeFaceArea(const OpenMesh::HalfedgeHandle heh);
	float computeVertexArea(const OpenMesh::VertexHandle vh);
	float cotan(double angle);
	std::pair<float, float> getOppositeAngles(const OpenMesh::HalfedgeHandle oh);
	float computeAngle(const OpenMesh::HalfedgeHandle heh);
	OpenMesh::Vec3f computeCentroid(const OpenMesh::HalfedgeHandle heh);
	Mesh::Normal filterFaceNormal(const OpenMesh::FaceHandle fh);

	float computeDistanceWeight(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle neighbour);

	OpenMesh::Vec3f cotangentLaplacian(const OpenMesh::VertexHandle vh);
	OpenMesh::Vec3f uniformLaplacian(const OpenMesh::VertexHandle vh);
	OpenMesh::Vec3f anisotropicLaplacian(const OpenMesh::VertexHandle vh);

	void taubinSmoothing(float lambda, float mu);
	void smoothMesh(int iterations);

	OpenMesh::Vec3f& laplacian_displacement(Mesh::VertexHandle _vh)
	{
		return mesh_.property(laplacian, _vh);
	}

private:
	OpenMesh::VPropHandleT<OpenMesh::Vec3f> laplacian;
	float variance_edge_length;

	float computeAverageEdgeLength();
	float computeVarianceEdgeLength();
};

