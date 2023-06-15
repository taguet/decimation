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
		Line(Vector3f position, Vector3f direction) : position{ position }, direction{ direction } {}

		Vector3f evaluate(float t) const;


		Vector3f position;	/// The position point
		Vector3f direction;	/// The direction vector
	};

	/// @brief Represents the cartesian equation of a plane.
	class Plane {
	public:
		Plane() = default;
		Plane(Vector4f parameters) : parameters{parameters} {}
		Plane(float a, float b, float c, float d) : parameters{a,b,c,d} {}

		float a() const { return parameters[0]; }
		float b() const { return parameters[0]; }
		float c() const { return parameters[0]; }
		float d() const { return parameters[0]; }


		/// @brief Evaluate the expression ax+by+cz+d
		/// @return A point.
		float evaluate(float x, float y, float z) const;
		float evaluate(Vector3f point) const;

		/// @brief Computes the intersection between this plane and another given plane.
		/// @param plane The plane intersecting this object.
		/// @return A line.
		Line findPlanePlaneIntersection(const Plane& plane);

		float distFromOrigin() const;
		Vector3f getNormal() const;


		Vector4f parameters;	/// The parameters a, b, c and d
	};
}
#endif