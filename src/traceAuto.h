#ifndef _TRACE_AUTO_H_
#define _TRACE_AUTO_H_

#include "typedefs.h"
std::vector<MyPolyline> traceAll(const cv::Mat & bwImg, const cv::Mat & origMask, const cv::Mat & extMask, const std::array<Eigen::MatrixXcd, 2>& roots, const Eigen::VectorXcd & X, const Eigen::MatrixXi& indices, std::map<std::array<int, 2>, std::vector<PixelInfo>>& pixelInfo, std::vector<std::array<bool, 2>>& endedWithASingularity);

#endif
