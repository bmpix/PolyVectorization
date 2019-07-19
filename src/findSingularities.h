#ifndef _FIND_SINGULARITIES_H_
#define _FIND_SINGULARITIES_H_
#include "typedefs.h"
#include <set>
bool isItASingularity(const Eigen::Vector2d & p, const std::array<std::complex<double>, 2> myRoots, const std::array<Eigen::MatrixXcd, 2>& roots, const cv::Mat & origMask);

std::set<std::array<int, 2>>  findSingularities(std::array<Eigen::MatrixXcd,2> & roots, Eigen::VectorXcd& X, const Eigen::MatrixXi& indices, const cv::Mat & origMask);

#endif
