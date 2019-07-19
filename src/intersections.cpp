#include "stdafx.h"
#include "intersections.h"

bool get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double * i_x, double * i_y, double * sOut, double * tOut)
{
	double s1_x, s1_y, s2_x, s2_y;
	s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
	s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

	double s, t;
	double sTop = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y));
	double sBottom = (-s2_x * s1_y + s1_x * s2_y);

	if ((sTop < 0 && sBottom > 0) || (sTop > 0 && sBottom < 0))
		return false;	

	if (sTop < 0)
	{
		sTop = -sTop;
		sBottom = -sBottom;
	}

	double tTop = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x));
	double tBottom = (-s2_x * s1_y + s1_x * s2_y);

	if ((tTop < 0 && tBottom > 0) || (tTop > 0 && tBottom < 0))
		return false;

	if (tTop < 0)
	{
		tTop = -tTop;
		tBottom = -tBottom;
	}

	if (sTop < sBottom && tTop < tBottom)
	{
		s =  sTop / sBottom;
		t =  tTop / tBottom;
		if (sOut != nullptr)
			*sOut = s;
		if (tOut != nullptr)
			*tOut = t;

		// Collision detected
		if (i_x != nullptr)
			*i_x = p0_x + (t * s1_x);
		if (i_y != nullptr)
			*i_y = p0_y + (t * s1_y);
		return true;
	}

	return false; // No collision
}
