#include "stdafx.h"
#include "polynomial_energy.h"
#include <numeric>
#include <iostream>
std::pair<double,Eigen::VectorXcd> polynomial_energy(const Eigen::VectorXcd & X, const Eigen::MatrixXd & weight, const Eigen::MatrixXcd & tau, const cv::Mat & mask, Eigen::MatrixXi& indices,std::vector<double>& energiesOut, bool computeGrad)
{
	int nonzeros = countNonZero(mask);
	Eigen::VectorXcd grad;
	if (computeGrad)
	{
		grad = Eigen::VectorXcd(nonzeros * 2);
		grad.setZero();
	}

	std::vector<double> energies(X.size()/2,0);
	double regEnergy = 0;
	int m = mask.rows;
	int n = mask.cols;
	
//#pragma omp parallel for
	for (int j=0; j<n; ++j)
		for (int i = 0; i<m; ++i)
		{
			double tmp = std::abs(tau(i, j));
			if ((mask.at<uchar>(i, j) == 0) || (std::abs(tau(i, j)) < 1e-6))
					continue;
			
			auto myTau = tau(i, j);
			int idx = indices(i, j);

			auto x0 = X(idx);
			auto x2 = X(idx + nonzeros);

			auto res = x0 + x2*std::pow(myTau,2) + std::pow(myTau,4);
			energies[idx] += (res.real()*res.real() + res.imag()*res.imag())*weight(i,j);

			if (computeGrad)
			{
				double a1 = x0.real(), b1 = x0.imag(), a2 = x2.real(), b2 = x2.imag();
				double r = myTau.real(), m = myTau.imag();
				
				grad(idx) = { 2 * a1 + (-2)*a2*pow(m,2) + 2 * pow(m,4) + (-4)*b2*m*r + 2 * a2*pow(r,2) + (-12)*pow(m,2)*pow(r,2) + 2 * pow(r,4),
					2 * b1 + (-2)*b2*pow(m,2) + 4 * a2*m*r + (-8)*pow(m,3)*r + 2 * b2*pow(r,2) + 8 * m*pow(r,3) };
				grad(idx) *= weight(i, j);

				grad(idx + nonzeros) = { (-2)*a1*pow(m,2) + 2 * a2*pow(m,4) + (-2)*pow(m,6) + 4 * b1*m*r + 2 * a1*pow(r,2) + 4 * a2*pow(m,2)*pow(r,2) + (-2)*pow(m,4)*pow(r,2) + 2 * a2*pow(r,4) + 2 * pow(m,2)*pow(r,4) + 2 * pow(r,6),
					(-2)*b1*pow(m,2) + 2 * b2*pow(m,4) + (-4)*a1*m*r + 4 * pow(m,5)*r + 2 * b1*pow(r,2) +
					4 * b2*pow(m,2)*pow(r,2) + 8 * pow(m,3)*pow(r,3) + 2 * b2*pow(r,4) + 4 * m*pow(r,5) };
				grad(idx+nonzeros) *= weight(i, j);
			}
		}

	double energy = std::accumulate(energies.begin(), energies.end(), 0.0);
	energiesOut = energies; //todo: get rid of the redundant copy
	return std::make_pair(energy, grad);
}
