#ifndef GEOMEQ_HPP
#define GEOMEQ_HPP

#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>

namespace Equation {
	using Vector3f = Eigen::Vector3f;
	using Vector4f = Eigen::Vector4f;

	typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;

	/// @brief Represents the parametric equation of a line.
	class Line {
	public:
		Line() = default;
		Line(const Vector3f& position, const Vector3f& direction) : origin{ position }, direction{ direction } {}

		Vector3f evaluate(const float t) const;

		Vector3f projectPoint(const Vector3f& p) const;
		Vector3f projectPoint(const Mesh::Point& p) const;
		float distToPoint(const Vector3f& p) const;
		float distToPoint(const Mesh::Point& p) const;


		Vector3f origin;	/// The position point
		Vector3f direction;	/// The direction vector
	};

	/// @brief Represents the cartesian equation of a plane.
	class Plane {
	public:
		Plane() = default;
		Plane(const Vector4f& parameters) : parameters{parameters} {}
		Plane(const float a, const float b, const float c, const float d) : parameters{a,b,c,d} {}

		float a() const { return parameters[0]; }
		float b() const { return parameters[1]; }
		float c() const { return parameters[2]; }
		float d() const { return parameters[3]; }


		/// @brief Evaluate the expression ax+by+cz+d
		/// @return A point.
		float evaluate(const float x, const float y, const float z) const;
		float evaluate(const Vector3f& point) const;
		float evaluate(const Mesh::Point& point) const;

		float solveX(const float y, const float z) const;
		float solveY(const float x, const float z) const;
		float solveZ(const float x, const float y) const;

		/// @brief Computes the intersection between this plane and another given plane.
		/// @param plane The plane intersecting this object.
		/// @return A line.
		Line findPlanePlaneIntersection(const Plane& plane) const;

		float distFromOrigin() const;
		float distToPoint(const Vector3f& p) const;
		float distToPoint(const Mesh::Point& p) const;
		float signedDistToPoint(const Vector3f& p) const;
		float signedDistToPoint(const Mesh::Point& p) const;
		Vector3f getNormal() const;
		Vector3f projectPoint(const Vector3f& p) const;
		Vector3f projectPoint(const Mesh::Point& p) const;

		bool isCoplanar(const Plane& plane) const;

		Vector4f parameters;	/// The parameters a, b, c and d
	};


	std::ostream& operator<<(std::ostream& os, const Plane& plane);
}
#endif