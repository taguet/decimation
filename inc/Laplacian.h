#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "MeshUtils.h"
#include "BilateralFiltering.h"


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
	
protected:
	BilateralFiltering filter;
};