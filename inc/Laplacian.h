#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "MeshUtils.h"


class Laplacian
{
public:
	using Mesh = OpenMesh::TriMesh_ArrayKernelT<>;

	Laplacian(Mesh& mesh);

	/// Compute the laplacian at a given vertex
	virtual OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) = 0;
	virtual void computeLaplacians();

	OpenMesh::Vec3f& laplacian_displacement(Mesh::VertexHandle _vh);

protected:
	Mesh* mesh_{ nullptr };
	OpenMesh::VPropHandleT<OpenMesh::Vec3f> laplacian;
};


class UniformLaplacian : public Laplacian
{
public:
	UniformLaplacian(Mesh& mesh);
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) override;
};


class CotangentLaplacian : public Laplacian
{
public:
	CotangentLaplacian(Mesh& mesh);
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) override;
};


class AnisotropicLaplacian : public Laplacian
{
public:
	AnisotropicLaplacian(Mesh& mesh);
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) override;
	void computeLaplacians();

	Mesh::Normal filterFaceNormal(const OpenMesh::FaceHandle fh, float threshold = 0.349066f);
	float computeDistanceWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2);
	float computeProximityWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2, float threshold = 0.349066f);

protected:
	Mesh::Normal weighFaceNormal(const OpenMesh::FaceHandle fh_i, const OpenMesh::FaceHandle fh_j, const float area, const float threshold);	

	float variance_edge_length;
};