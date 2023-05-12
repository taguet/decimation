#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "MeshUtils.h"


class Laplacian
{
public:
	using Mesh = OpenMesh::TriMesh_ArrayKernelT<>;

	Laplacian(Mesh& mesh);

	/// @brief Computes the laplacian displacement at a given vertex.
	/// 
	/// This function is implemented by the derived classes.
	/// 
	/// @param vh the handle of the vertex
	/// @return A vector of floats corresponding to the computed displacement for that vertex.
	virtual OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) = 0;

	/// @brief Computes the laplacian displacement of all vertices and stores them in a property.
	virtual void computeLaplacians();

	/// @brief Getter for the laplacian displacement
	/// @param _vh the given vertex handle
	/// @return A vector stored in the laplacian displacement property.
	OpenMesh::Vec3f& laplacian_displacement(Mesh::VertexHandle _vh);

protected:
	Mesh* mesh_{ nullptr };	///< Pointer to the mesh to do all computations on.
	OpenMesh::VPropHandleT<OpenMesh::Vec3f> laplacian;	///< Property storing the laplacian displacement on each vertex.
};


/// @brief An implementation of the uniform laplacian operator.
class UniformLaplacian : public Laplacian
{
public:
	UniformLaplacian(Mesh& mesh);
	/// @brief Computes the laplacian displacement at a given vertex.
	/// 
	/// This displacement is given by \f[
	///		\Delta _{v_i}=\frac{1}{\left | N_{v_i} \right |} \sum_{v_j\in N_{v_i}}(v_j-v_i)
	/// \f]
	/// 
	/// @param vh A given vertex handle.
	/// @return A vector of floats corresponding to the computed displacement for that vertex.
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) override;
};


/// @brief An implementation of the cotangent approximation of the Laplace-Beltrami operator.
class CotangentLaplacian : public Laplacian
{
public:
	CotangentLaplacian(Mesh& mesh);
	/// @brief Computes the laplacian displacement at a given vertex.
	/// 
	/// This displacement is given by \f[
	///		\Delta _{v_i}=\frac{1}{2A_i} \sum_{v_j\in N_{v_i}} [cot(\alpha_{ij}) + cot(\beta_{ij})](v_j-v_i)
	/// \f]
	/// 
	/// @param vh A given vertex handle.
	/// @return A vector of floats corresponding to the computed displacement for that vertex.
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) override;
};


/// @brief An implementation of the discrete anisotropic laplacian.
class AnisotropicLaplacian : public Laplacian
{
public:
	AnisotropicLaplacian(Mesh& mesh);
	/// @brief Computes the laplacian displacement at a given vertex.
	/// 
	/// This displacement is given by \f[
	///		\Delta _{v_i}=\frac{1}{\left | N_{v_i}(f)\right |_{0}} 
	///		\sum_{f_{k} \in N_{v_{i}}(f)} [(c_{f_{k}}-v_i)\cdot {\hat n}'_{f_{k}}]\cdot n_{f_{k}} 
	/// \f]
	/// 
	/// @param vh A given vertex handle.
	/// @return A vector of floats corresponding to the computed displacement for that vertex.
	OpenMesh::Vec3f computeLaplacian(const OpenMesh::VertexHandle vh) override;
	void computeLaplacians();

	/// @brief Computes the filtered normal of a given face.
	/// @param fh A given face handle.
	/// @param threshold The maximum angle between two normals before considering them as not belonging to the same neighbourhood.
	/// @return A vector representing the filtered normal.
	Mesh::Normal filterFaceNormal(const OpenMesh::FaceHandle fh, float threshold = 0.349066f);

	/// @brief Computes the spatial distance-based weight between two faces
	/// 
	/// This weight is given by \f[
	///		\alpha_{f_if_j}=exp\left (-\frac{\left \| c_{f_i}-c_{f_j} \right \|^2}{2\sigma_{dist}^2} \right )
	/// \f],
	/// with \f[\sigma_{dist}\f] the variance parameter of distance proximity, estiated by the average edge length of the mesh edges.
	/// 
	/// @param fh1 A given face handle. 
	/// @param fh2 The neighbouring face's handle.
	/// @return A scalar between 0 and 1. The closer those faces are, the closer the scalar is to 1.
	float computeDistanceWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2);

	/// @brief Computes the spatial distance-based weight between two faces
	/// 
	/// This weight is given by \f[
	///		\beta_{f_if_j}=exp\left (-\frac{\left \| 1-n_{f_i}\cdot n_{f_j} \right \|^2}{(1-cos(\theta))^2} \right )
	/// \f],
	/// with \f[\theta\f] the angle threshold, set by default to 20 degrees.
	/// 
	/// @param fh1 A given face handle. 
	/// @param fh2 The neighbouring face's handle.
	/// @return A scalar between 0 and 1. The closer more similar those normals are, the closer the scalar is to 1.
	float computeProximityWeight(const OpenMesh::FaceHandle fh1, const OpenMesh::FaceHandle fh2, float threshold = 0.349066f);

protected:
	Mesh::Normal weighFaceNormal(const OpenMesh::FaceHandle fh_i, const OpenMesh::FaceHandle fh_j, const float area, const float threshold);	

	float variance_edge_length;
};