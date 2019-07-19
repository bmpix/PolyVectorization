#ifndef _REMOVE_SHORT_BRANCHES_H_
#define _REMOVE_SHORT_BRANCHES_H_

#include "typedefs.h"
#include "graph_typedefs.h"

void removeBranchesFilter1(G& g, bool onlyAtIntersections, const std::map<edge_descriptor, size_t>& myChainsBeforeSharpCornersSplit);

#endif