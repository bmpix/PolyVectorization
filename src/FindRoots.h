#ifndef _FIND_ROOTS_H_
#define _FIND_ROOTS_H_
#include <array>
#include "Eigen/Dense"
#include "opencv2/core/core.hpp"

std::array<std::complex<double>, 2> findRoots(std::complex<double> c0, std::complex<double> c2);
std::array<Eigen::MatrixXcd, 2> findRoots(const Eigen::VectorXcd& X, const cv::Mat& mask);

#endif
