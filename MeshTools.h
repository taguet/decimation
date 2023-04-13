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
	double computeFaceArea(const OpenMesh::HalfedgeHandle heh);
	double computeVertexArea(const OpenMesh::VertexHandle vh);
	double cotan(double angle);
	std::pair<double, double> getOppositeAngles(const OpenMesh::HalfedgeHandle oh);

	template <class T>
	T discreteLaplacian(const OpenMesh::VertexHandle vh, T (*function)(const OpenMesh::VertexHandle));
};

