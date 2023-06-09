#pragma once

#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <list>
#include <set>
#include "GeomEq.hpp"

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;
using Vector3f = Eigen::Vector3f;
using Vector4f = Eigen::Vector4f;
using VectorXf = Eigen::VectorXf;
using MatrixXf = Eigen::MatrixXf;
using Matrix2f = Eigen::Matrix2f;

/// @brief Utility class with basic operations to use in algorithms on a mesh.
class MeshUtils
{
public:
	/// @brief Computes the area of a given face.
	/// @param mesh
	/// @param heh A halfedge handle belonging to the face.
	/// @return The area.
	static float computeFaceArea(Mesh& mesh, const OpenMesh::HalfedgeHandle heh);
	static float computeFaceArea(Mesh& mesh, const OpenMesh::FaceHandle fh);

	/// @brief Computes the vertex area of a given vertex.
	/// 
	/// The vertex area is one third of the sum of the areas of all faces adjacent to the vertex.
	/// 
	/// @param mesh 
	/// @param vh The vertex handle
	/// @return The vertex area.
	static float computeVertexArea(Mesh& mesh, const OpenMesh::VertexHandle vh);


	/// @brief Computes the opposite angles on both sides of an edge.
	/// @param mesh 
	/// @param oh The halfedge handle
	/// @return Two angles in radians.
	static std::pair<float, float> getOppositeAngles(Mesh& mesh, const OpenMesh::HalfedgeHandle oh);

	/// @brief Computes the angle between a halfedge and its successor.
	/// @param mesh 
	/// @param heh The halfedge handle.
	/// @return An angle in radians.
	static float computeAngle(Mesh& mesh, const OpenMesh::HalfedgeHandle heh);

	static float cotan(float angle);

	/// @brief Computes the centroid of a face.
	/// @param mesh 
	/// @param heh An halfedge handle belonging to the face.
	/// @return A vector of floats representing the position of the centroid.
	static OpenMesh::Vec3f computeCentroid(Mesh& mesh, const OpenMesh::HalfedgeHandle heh);


	/// @brief Computes the average of the lengths of all edges.
	/// @param mesh 
	/// @return The average length.
	static float computeAverageEdgeLength(Mesh& mesh);

	/// @brief Computes the variance of the lengths of all edges.
	/// @param mesh 
	/// @return The variance of the edges' length.
	static float computeVarianceEdgeLength(Mesh& mesh);


	static void getFaceNeighbors(Mesh& mesh, const Mesh::FaceHandle, std::list<Mesh::FaceHandle>& neighbors);


	/// @brief Projected a given vertex onto a line.
	/// @param mesh 
	/// @param vh 
	/// @param line A parametric equation of a line.
	static void projectVertexToLine(Mesh& mesh, const Mesh::VertexHandle vh, const Equation::Line& line);


	/// @brief Fits plane to a set of vertices by minimizing the orthogonal distance by least squares.
	/// @param mesh 
	/// @param vertices 
	/// @return The cartesian equation of the fitted plane.
	static Equation::Plane fitPlaneToVertices(Mesh& mesh, std::set<Mesh::VertexHandle>& vertices); 
};

