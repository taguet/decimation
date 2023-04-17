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
	double computeFaceArea(const OpenMesh::HalfedgeHandle heh);
	double computeVertexArea(const OpenMesh::VertexHandle vh);
	float cotan(double angle);
	std::pair<float, float> getOppositeAngles(const OpenMesh::HalfedgeHandle oh);

	OpenMesh::Vec3f discreteLaplacian(const OpenMesh::VertexHandle vh);
	OpenMesh::Vec3f uniformLaplacian(const OpenMesh::VertexHandle vh);

	void smoothMesh(int h, float lambda);
	void smoothMesh(int iterations);

	OpenMesh::Vec3f& laplacian_displacement(Mesh::VertexHandle _vh)
	{
		return mesh_.property(laplacian, _vh);
	}

private:
	OpenMesh::VPropHandleT<OpenMesh::Vec3f> laplacian;
};

