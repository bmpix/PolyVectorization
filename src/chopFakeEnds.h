#ifndef _CHOP_FAKE_ENDS_H_
#define _CHOP_FAKE_ENDS_H_

#include "typedefs.h"
#include "graph_typedefs.h"

std::pair<std::vector<MyPolyline>,G> chopFakeEnds(const std::vector<MyPolyline>& polys, const std::vector<std::vector<double>>& radii, const std::vector<std::array<bool, 2>>& protectedEnds, 
	const std::vector<std::array<bool, 2>>& isItASpecialDeg2Vertex, const std::vector<std::pair<PointOnCurve, PointOnCurve>>& yJunctions);

#endif
