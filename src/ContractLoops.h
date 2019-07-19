#ifndef _CONTRACT_LOOPS_H_

#include "typedefs.h"
#include "graph_typedefs.h"

std::vector<edge_descriptor> contractLoops(G& reebGraph, const cv::Mat & origMask, const std::vector<MyPolyline>& polys); //returns incontractible loops


#endif