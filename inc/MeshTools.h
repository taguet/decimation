#pragma once

#include "gl.h"
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "Laplacian.h"


typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class MeshTools
{
private:
	Mesh* mesh_{ nullptr };
public:
	
	MeshTools(void);
	MeshTools(Mesh& mesh);
	void setMesh(Mesh &mesh);

	//void taubinSmoothing(float lambda, float mu, int iterations=1);
	void smoothMesh(Laplacian& laplacian, int iterations=1, float factor=0.2);
};

