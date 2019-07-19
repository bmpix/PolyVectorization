#ifndef _TOTAL_ENERGY_H_
#define _TOTAL_ENERGY_H_
#include "Eigen/Dense"
#include "opencv2/core/core.hpp"
#include <map>
#include <array>
struct TotalEnergy
{
	TotalEnergy(cv::Mat & bwImg, const Eigen::MatrixXd & weight, const Eigen::MatrixXcd & tauNormalized, double beta, cv::Mat & mask, const Eigen::MatrixXi& indices, int nnz);
	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
	double energyOnly(const Eigen::VectorXd & x);

	int m, n;
	const double beta;
	double alpha;
	cv::Mat& bwImg;
	const Eigen::MatrixXd& weight;
	const Eigen::MatrixXcd & tauNormalized;
	Eigen::MatrixXcd g;
	cv::Mat & mask;
	bool noisy;
	Eigen::MatrixXd smartWeights;
	Eigen::MatrixXi indices;
	int nnz;
	Eigen::MatrixXd onesMatrix;
};

#endif