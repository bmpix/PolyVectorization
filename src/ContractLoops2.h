#ifndef _CONTRACT_LOOPS2_H_

#include "typedefs.h"
#include "graph_typedefs.h"
#include <stack>

std::vector<edge_descriptor> contractLoops2(G& g, const cv::Mat & origMask, const std::vector<MyPolyline>& polys); //returns incontractible loops


#endif
