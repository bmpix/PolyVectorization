#ifndef _L2_REGULARIZER_H_
#define _L2_REGULARIZER_H_
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "opencv2/core/core.hpp"

std::pair<double, Eigen::VectorXcd> l2_regularizer(const Eigen::VectorXcd& X, const Eigen::MatrixXd& D, const cv::Mat & mask, Eigen::MatrixXi& indices, bool computeGrad);
Eigen::SparseMatrix<double> laplacian_matrix(const cv::Mat& mask, const Eigen::MatrixXi& indices, const Eigen::VectorXd& w);

#endif