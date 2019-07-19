#ifndef _IS_LOOP_CONTRACTIBLE_H_
#define _IS_LOOP_CONTRACTIBLE_H_

#include "typedefs.h"
#include "graph_typedefs.h"

bool isLoopContractible(const std::vector<edge_descriptor>& loop, const cv::Mat& origMask, const G& g, const std::vector<MyPolyline>& polys, std::vector<cv::Point2f>& outCvLoop);

#endif