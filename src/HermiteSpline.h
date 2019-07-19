#ifndef _HERMITE_SPLINE_H_
#define _HERMITE_SPLINE_H_

#include "typedefs.h"
#include <array>

class HermiteSpline
{
public:
	HermiteSpline::HermiteSpline(const Eigen::Vector2d& P0, const Eigen::Vector2d& N0, const Eigen::Vector2d& P1, const Eigen::Vector2d& N1)
		:P0_(P0), P1_(P1), N0_(N0), N1_(N1)
	{
	}

	/* Centripetal Catmull-Rom: interpolating spline for four points. Draws from P1 to P2 only.
	Reference: http://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
	*/
	static HermiteSpline fromCatmullRom(const std::array<Eigen::Vector2d, 4>& P)
	{
		/*Taken from http://stackoverflow.com/a/23980479/2929337 */

		std::array<double, 4> t;
		std::array<double, 3> dt;

		for (int i = 0; i < 3; i++)
		{
			dt[i] = std::pow((P[i + 1] - P[i]).squaredNorm(), 0.25);
		}

		// safety check for repeated points
		if (dt[1] < 1e-7)    dt[1] = 1.0;
		if (dt[0] < 1e-7)    dt[0] = dt[1];
		if (dt[2] < 1e-7)    dt[2] = dt[1];

		t[0] = 0;
		for (int i = 1; i < 4; i++)
			t[i] = t[i - 1] + dt[i - 1];

		Eigen::Vector2d N1 = (P[1] - P[0]) / dt[0] - (P[2] - P[0]) / (dt[0] + dt[1]) + (P[2] - P[1]) / dt[1];
		Eigen::Vector2d N2 = (P[2] - P[1]) / dt[1] - (P[3] - P[1]) / (dt[1] + dt[2]) + (P[3] - P[2]) / dt[2];
		N1 *= dt[1];
		N2 *= dt[2];
		return HermiteSpline(P[1], N1, P[2], N2);
	}

	double HermiteSpline::getLength() const
	{
		//The length of Hermite Curve is can't be calculated analytically.
		//Instead of this, one may use Runge-Cutta method, but this works as well.
		Eigen::Vector2d prevPoint = getPoint(0);
		const int samples = 50;
		double length = 0;
		for (int i = 1; i <= samples; i++)
		{
			double t = (double)i / samples;
			Eigen::Vector2d p = getPoint(t);
			length += (p - prevPoint).norm();
			prevPoint = p;
		}
		return length;
	}


	Eigen::Vector2d HermiteSpline::getPoint(double t) const
	{
		return r(t);
	}


	Eigen::Vector2d HermiteSpline::getTangent(double t) const
	{
		return dr(t);
	}

	double HermiteSpline::getAbsCurvature(double t)
	{
		Eigen::Vector2d drt = dr(t);
		Eigen::Vector2d ddrt = ddr(t);
		double drAbs = drt.norm();
		Eigen::Vector3d dr3(drt.x(), drt.y(), 0);
		Eigen::Vector3d ddr3(ddrt.x(), ddrt.y(), 0);
		return dr3.cross(ddr3).norm() / (drAbs*drAbs*drAbs);
	}

	double HermiteSpline::integrateAbsCurvature()
	{
		const int nSteps = 100;
		double sum = 0;
		double max = 0;
		double arcLength = 0;
		Eigen::Vector2d p = r(0);
		for (int i = 0; i < nSteps; i++)
		{
			double t = (double)i / nSteps;
			double k = getAbsCurvature(t);
			Eigen::Vector2d newP = r(t);
			sum += k*(newP - p).norm();
			p = newP;
		}

		return sum;
	}



	HermiteSpline& operator =(const HermiteSpline &p)
	{
		P0_ = p.P0_;
		N0_ = p.N0_;
		P1_ = p.P1_;
		N1_ = p.N1_;
		return *this;
	}

private:
	//radius-vector and its derivatives
	Eigen::Vector2d HermiteSpline::r(double t) const
	{
		double t3 = t*t*t;
		double t2 = t*t;
		return P0_*(2 * t3 - 3 * t2 + 1) + N0_*(t3 - 2 * t2 + t) + P1_*(-2 * t3 + 3 * t2) + N1_*(t3 - t2);
	}

	Eigen::Vector2d HermiteSpline::dr(double t) const
	{
		double t2 = t*t;
		return P0_*(6 * t2 - 6 * t) + N0_*(3 * t2 - 4 * t + 1) + P1_*(-6 * t2 + 6 * t) + N1_*(3 * t2 - 2 * t);
	}

	Eigen::Vector2d HermiteSpline::ddr(double t) const
	{
		return P0_*(12 * t - 6) + N0_*(6 * t - 4) + P1_*(-12 * t + 6) + N1_*(6 * t - 2);
	}

private:
	Eigen::Vector2d P0_, P1_;
	Eigen::Vector2d N0_, N1_;
};

#endif