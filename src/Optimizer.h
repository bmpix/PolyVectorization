#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "opencv2/core/core.hpp"

Eigen::VectorXcd optimize(cv::Mat& bwImg, const Eigen::MatrixXd& weight, const Eigen::MatrixXcd& tauNormalized, double beta, cv::Mat& extMask, const Eigen::MatrixXi& indices);

#endif
