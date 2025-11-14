#ifndef _POLYNOMIAL_ENERGY_H_
#define _POLYNOMIAL_ENERGY_H_
#include "Eigen/Dense"
#include "opencv2/core/core.hpp"
#include "Eigen/Sparse"

std::pair<double, Eigen::VectorXcd> polynomial_energy(const Eigen::VectorXcd& X, const Eigen::MatrixXd & weight, const Eigen::MatrixXcd& tau, const cv::Mat & mask, Eigen::MatrixXi& indices, std::vector<double>& energiesOut, bool computeGrad);
std::pair<Eigen::SparseMatrix<std::complex<double>>, Eigen::VectorXcd> polynomial_energy_matrix(const Eigen::MatrixXd& weight, const Eigen::MatrixXcd& tau, const cv::Mat& mask, const Eigen::MatrixXi& indices);

#endif
