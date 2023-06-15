#ifndef GEOMEQ_HPP
#define GEOMEQ_HPP

#include <Eigen/Dense>

namespace Equation {
	using Vector3f = Eigen::Vector3f;
	using Vector4f = Eigen::Vector4f;

	/// @brief Represents the parametric equation of a line
	class Line {
	public:
		Line(Vector3f position, Vector3f direction) : position{ position }, direction{ direction } {}

		Vector3f position;	/// The position point
		Vector3f direction;	/// The direction vector
	};

	/// @brief Represents the cartesian equation of a plane
	class Plane {
	public:
		Plane(Vector3f parameters) : parameters{parameters} {}
		Plane(float a, float b, float c, float d) : parameters{a,b,c,d} {}

		Vector4f parameters;	/// The parameters a, b, c and d
	};
}
#endif