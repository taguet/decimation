#include "GeomEq.hpp"
#include <iostream>

namespace Equation {
	//------------------- PLANE --------------------//

	/// Uses an approach from Gellert et al. 1989, p.542 to find a parametric equation of the line.
	Line Plane::findPlanePlaneIntersection(const Plane& plane) const {
		Vector3f n1{ getNormal() };
		Vector3f n2{ plane.getNormal() };
		Eigen::MatrixXf m{ 3,2 };	// Matrix m = [n1 n2]^T
		m.col(0) = n1;
		m.col(1) = n2;
		m.transposeInPlace();
		Eigen::VectorXf b{ 2 };		// Vector b = -[p1 p2]^T
		b << -distFromOrigin(), -plane.distFromOrigin();

		Vector3f point{ m.colPivHouseholderQr().solve(b)};	// mx = b
		Eigen::FullPivLU<Eigen::MatrixXf> lu{ m };
		Eigen::MatrixXf m_null_space{ lu.kernel() };
		Vector3f direction{ m_null_space.col(0)};	//mx = 0 to find a vector given by the nullspace of m
		return { point, direction };
	}


	Vector3f Plane::getNormal() const {
		Vector3f normal{ a(), b(), c() };
		return normal.normalized();
	}


	Vector3f Plane::projectPoint(const Vector3f& p) const {
		float dist{ signedDistToPoint(p) };
		Vector3f normal{ getNormal() };
		return p - dist * normal;
	}


	Vector3f Plane::projectPoint(const Mesh::Point& p) const {
		return projectPoint(Vector3f{ p[0], p[1], p[2] });
	}


	bool Plane::isCoplanar(const Plane& plane) const {
		Vector3f n0{ getNormal() };
		Vector3f n1{ plane.getNormal() };
		std::cerr << n0.dot(n1) << '\n';
		return n0.dot(n1) >= cos(0.349066f);
	}


	float Plane::distFromOrigin() const {
		Vector3f normal{ getNormal() };
		return d() / normal.norm();
	}


	float Plane::distToPoint(const Vector3f& p) const {
		return std::abs(signedDistToPoint(p));
	}


	float Plane::distToPoint(const Mesh::Point& p) const {
		return distToPoint(Vector3f{p[0], p[1], p[2]});
	}


	float Plane::signedDistToPoint(const Vector3f& p) const {
		return (a() * p[0] + b() * p[1] + c() * p[2] + d()) / std::sqrtf(a() * a() + b() * b() + c() * c());
	}


	float Plane::signedDistToPoint(const Mesh::Point& p) const {
		return signedDistToPoint(Vector3f{ p[0], p[1], p[2] });
	}


	float Plane::evaluate(const float x, const float y, const float z) const {
		Eigen::Vector4f X{ x,y,z,1 };
		auto result{ parameters.transpose() * X };
		return result.eval()(0,0);
	}
	
	
	float Plane::evaluate(const Vector3f& point) const {
		Eigen::Vector4f X{ point(0), point(1), point(2), 1};
		auto result{ parameters.transpose() * X };
		return result.eval()(0,0);
	}


	float Plane::evaluate(const Mesh::Point& point) const {
		return evaluate(point[0], point[1], point[2]);
	}


	float Plane::solveX(const float y, const float z) const {
		return -(b() * y + c() * z + d()) / a();
	}


	float Plane::solveY(const float x, const float z) const {
		return -(a() * x + c() * z + d()) / b();
	}


	float Plane::solveZ(const float x, const float y) const {
		return -(a() * x + b() * y + d()) / c();
	}

	std::ostream& operator<<(std::ostream& os, const Plane& plane) {
		os << plane.a() << "x + " << plane.b() << "y + " << plane.c() << "z + " << plane.d() << " = 0";
		return os;
	}

//------------------- LINE --------------------//

	Vector3f Line::evaluate(const float t) const {
		return origin + t * direction;
	}


	Vector3f Line::projectPoint(const Vector3f& p) const {
		return (direction.dot(p) / direction.dot(direction)) * direction + origin;
	}


	Vector3f Line::projectPoint(const Mesh::Point& p) const {
		return projectPoint(Vector3f{ p[0], p[1], p[2] });
	}


	float Line::distToPoint(const Vector3f& p) const {
		Vector3f proj{ projectPoint(p) };
		return p.norm();
	}


	float Line::distToPoint(const Mesh::Point& p) const {
		return distToPoint(Vector3f{ p[0], p[1], p[2] });
	}
}