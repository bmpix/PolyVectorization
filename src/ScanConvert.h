#pragma once
#include "graph_typedefs.h"

std::set<std::array<int, 2>> scanConvert(const std::vector<Eigen::Vector2d>& poly, const cv::Mat & origMask, const std::array<Eigen::MatrixXcd, 2>& roots,
	int startIdx, std::map<std::array<int, 2>, bool>& isRootMatchingOK, std::array<bool, 2>&  hitSingularity,
	const std::set<std::array<int, 2>>& singularities, std::array<double, 2> segment);