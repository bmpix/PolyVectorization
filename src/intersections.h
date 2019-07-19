#pragma once
#include <numeric>
#include <algorithm>

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the doubles i_x and i_y.
bool get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
	double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y, double *sOut, double *tOut);