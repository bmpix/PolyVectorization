#ifndef _GREEDY_TRACE_H_
#define _GREEDY_TRACE_H_

#include "typedefs.h"
#include <set>
std::pair<MyPolyline, std::array<bool, 2>> greedyTrace(const cv::Mat & origMask, 
	const std::array<Eigen::MatrixXcd, 2>& allRoots, const Eigen::Vector2d & seedCenter, const Eigen::Vector2d& initialDir, const Eigen::VectorXcd & X, 
	const std::map<std::array<int, 2>, std::vector<PixelInfo>>& pixelInfo, std::map<std::array<int, 2>, std::vector<PixelInfo>>& outNewPixelInfo, 
	const Eigen::MatrixXi& indices, int myCurveIdx, bool perpendicular,  std::array<double,2>& closestDistToSingularity);

#endif
