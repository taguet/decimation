#ifndef GEOMEQ_HPP
#define GEOMEQ_HPP

#include <Eigen/Dense>

namespace Equation {
	using Vector3f = Eigen::Vector3f;
	using Vector4f = Eigen::Vector4f;

	/// @brief Represents the parametric equation of a line.
	class Line {
	public:
		Line() = default;
		Line(const Vector3f& position, const Vector3f& direction) : position{ position }, direction{ direction } {}

		Vector3f evaluate(const float t) const;

		Vector3f projectPoint(const Vector3f& p) const;
		float distToPoint(const Vector3f& p) const;


		Vector3f position;	/// The position point
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
		float evaluate(const Vector3f point) const;

		/// @brief Computes the intersection between this plane and another given plane.
		/// @param plane The plane intersecting this object.
		/// @return A line.
		Line findPlanePlaneIntersection(const Plane& plane);

		float distFromOrigin() const;
		float distToPoint(const Vector3f& p) const;
		float signedDistToPoint(const Vector3f& p) const;
		Vector3f getNormal() const;
		Vector3f projectPoint(const Vector3f& p) const;


		Vector4f parameters;	/// The parameters a, b, c and d
	};
}
#endif