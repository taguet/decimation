#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class Laplacian
{
public:
	Laplacian(Mesh& mesh);
	~Laplacian();

	/// Compute the laplacian at a given vertex
	virtual OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) = 0;

protected:
	Mesh* mesh_{ nullptr };
};


class UniformLaplacian : public Laplacian
{
public:
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh);
};


class CotangentLaplacian : public Laplacian
{
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh);
};


class AnisotropicLaplacian : public Laplacian
{
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh);
};