#include "GeomEq.hpp"

namespace Equation {
	/// Uses an approach from Gellert et al. 1989, p.542 to find a parametric equation of the line.
	Line Plane::findPlanePlaneintersection(const Plane& plane) {
		Vector3f n1{ getNormal() };
		Vector3f n2{ plane.getNormal() };
		Eigen::MatrixXf m{ 3,3 };	// Matrix m = [n1 n2]^T
		m.col(0) = n1;
		m.col(1) = n2;
		m = m.transpose();
		Eigen::VectorXf b{ 2 };		// Vector b = -[p1 p2]^T
		b << -distFromOrigin(), -plane.distFromOrigin();

		Vector3f point{ m.colPivHouseholderQr().solve(b)};	// mx = b
		Vector3f direction{ m.colPivHouseholderQr().solve(Vector3f::Zero()) };	//mx = 0 to find a vector given by the nullspace of m
		return { point, direction };
	}


	Vector3f Plane::getNormal() const {
		Vector3f normal{ a(), b(), c() };
		return normal.normalized();
	}


	float Plane::distFromOrigin() const {
		Vector3f normal{ getNormal() };
		return d() / normal.norm();
	}
}