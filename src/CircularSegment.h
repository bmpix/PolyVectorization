#ifndef _CIRCULAR_SEGMENT_H_
#define _CIRCULAR_SEGMENT_H_

/*
 Defines operations for a segment on a (potentially) closed curve.
 Examples of valid segments are: [1.5, 5] (for any modulus), [9, 0.5] for modulus 10.
 Modulus = length of the curve or number of pieces in the polyline.
*/
#include <vector>

class CircularSegment 
{
public:
	//if modulus <0, the curve is considered to be NOT closed
	CircularSegment (double begin=-1, double end=-1, double modulus=-1);

	//if intersection is empty, returns FALSE
	bool intersectionWith (const CircularSegment& other,CircularSegment& outIntersection, double tolerance = 1e-10) const;

	//union is not always a segment. in case it's not a segment (i.e. two disjoint segments or empty set), the function returns FALSE
	std::vector<CircularSegment> unionWith (const CircularSegment& other, double tolerance = 1e-10) const;

	std::vector<CircularSegment> complement() const;

	void extendToContainPoint (double p);
	bool isInside (double p, double tolerance = 1e-10) const;
	double length() const;
	double pointAtT(double t) const;
	const double& begin() const {return begin_;}
	double& begin() {return begin_;}
	const double& end() const {return end_;}
	double& end() {return end_;}
	const double& modulus() const {return modulus_;}
	double& modulus() {return modulus_;}
private:
	double begin_, end_, modulus_;
	bool empty_;
};

#endif