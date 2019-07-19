#ifndef _BILINEAR_INTERP_COEFFS_H_
#define _BILINEAR_INTERP_COEFFS_H_
#include "Eigen/Dense"
#include "opencv2/core/core.hpp"
#include <array>

std::array<std::complex<double>, 2> bilinearInterpCoeff(const Eigen::VectorXcd &X, const Eigen::Vector2d& p, const cv::Mat& extMask, const Eigen::MatrixXi& indices);

#endif
