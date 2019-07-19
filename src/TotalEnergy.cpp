#include "stdafx.h"
#include "TotalEnergy.h"
#include "polynomial_energy.h"
#include "l2_regularizer.h"
#include <iostream>
#include <fstream>
#include "Params.h"
TotalEnergy::TotalEnergy(cv::Mat & bwImg, const Eigen::MatrixXd & weight, const Eigen::MatrixXcd & tauNormalized, double beta, cv::Mat & mask, const Eigen::MatrixXi& indices, int nnz)
	:bwImg(bwImg),weight(weight),tauNormalized(tauNormalized),beta(beta),mask(mask),noisy(noisy),indices(indices),nnz(nnz)
{
	m = bwImg.rows;
	n = bwImg.cols;

	alpha = FRAME_FIELD_REGULARIZER_WEIGHT;
	smartWeights = Eigen::MatrixXd(m,n);
	smartWeights.setOnes();
	for (int i=1; i<m-1; ++i)
		for (int j = 1; j < n-1; ++j)
		{
			smartWeights(i, j) = 0.25*((1-weight(i - 1, j)) + (1-weight(i + 1, j)) + (1-weight(i, j - 1)) + (1-weight(i, j + 1)));
		}

	onesMatrix = Eigen::MatrixXd(m, n);
	onesMatrix.setOnes();

	g = tauNormalized * std::complex<double>(0, 1);
}

double TotalEnergy::operator()(const Eigen::VectorXd & x, Eigen::VectorXd & grad)
{
	Eigen::VectorXcd g1, g2, g3;
	double e1, e2, e3;
	Eigen::VectorXcd x_complex = x.head(x.size() / 2) + std::complex<double>(0, 1)*x.tail(x.size() / 2);
	std::vector<double> energiesOut;
	std::tie(e1,g1) = polynomial_energy(x_complex, weight, tauNormalized, mask, indices, energiesOut,true);

	std::tie(e2, g2) = polynomial_energy(x_complex, onesMatrix, g, mask, indices, energiesOut, true);
	int c = cv::countNonZero(mask);
	std::tie(e3, g3) = l2_regularizer(x_complex, smartWeights, mask, indices,true);

	double totalEnergy = e1 + alpha*e2 + beta*e3;
	Eigen::VectorXcd grad_complex = g1 + alpha*g2 + beta*g3;
	grad.setZero();
	grad << grad_complex.real() , grad_complex.imag();

	return totalEnergy;
}

double TotalEnergy::energyOnly(const Eigen::VectorXd & x)
{
	Eigen::VectorXcd g1, g2, g3;
	double e1, e2, e3;
	Eigen::VectorXcd x_complex = x.head(x.size() / 2) + std::complex<double>(0, 1)*x.tail(x.size() / 2);
	std::vector<double> energiesOut;
	std::tie(e1, g1) = polynomial_energy(x_complex, weight, tauNormalized, mask, indices, energiesOut, false);

	std::tie(e2, g2) = polynomial_energy(x_complex, onesMatrix, g, mask, indices, energiesOut, false);
	std::tie(e3, g3) = l2_regularizer(x_complex, smartWeights, mask, indices, false);

	double totalEnergy = e1 + alpha*e2 + beta*e3;
	return totalEnergy;
}