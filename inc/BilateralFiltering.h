#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "MeshUtils.h"

class BilateralFiltering
{
public:
	BilateralFiltering(Mesh& mesh);
	~BilateralFiltering();

	void initializeNormals();

	void filterFaceNormals(int iteration = 1, float threshold = 0.349066f);

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

	void update();

	Mesh::Normal& filteredNormal(const Mesh::FaceHandle fh);

protected:
	Mesh::Normal weighFaceNormal(const OpenMesh::FaceHandle fh_i, const OpenMesh::FaceHandle fh_j, const float area, const float threshold);
	float variance_edge_length;
	Mesh* mesh{ nullptr };
	OpenMesh::FPropHandleT<Mesh::Normal> filtered_normals;
};

