#pragma once


#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <algorithm> 
#include <fstream>

#define TWOPI 6.2831853
#define PI 3.14159268
#define EPS 10e-4

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

class ComputingTools
{

private:

	float volume_; 
	float area_; 
	float mean_curv_glob_;
	float gauss_curv_glob_;
	Mesh mesh_;

	std::vector<std::vector<float>> principal_curvatures_;

public:

	ComputingTools(void);

	ComputingTools(Mesh _mesh);

	virtual ~ComputingTools(void);

	double volume();

	void compute_volume();

	double area();

	void compute_area();

	OpenMesh::Vec3f halfedgeVector(Mesh::HalfedgeHandle const &h) const;

	void compute_local_principalCurvatures(float &_k1, float &_k2, const OpenMesh::VertexHandle _vh);

	void set_mesh(Mesh& _mesh);

	std::vector<std::vector<float>> principal_curvatures();

	void calc_princ_curvatures();

	float cot(float _t);
};
