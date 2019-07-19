#include "stdafx.h"
#include "Simplify.h"
#include <queue>

typedef std::pair<Eigen::Vector2d, Eigen::Vector2d> Segment;

double perpendicularDistance(const Eigen::Vector2d& point, Segment segment)
{
	Eigen::Vector2d v = segment.first - segment.second;
	Eigen::Vector2d w = point - segment.second;

	Eigen::Vector2d closestPoint;

	double c1 = w.dot(v);

	if (c1<0)
	{
		closestPoint = segment.second;
		return (closestPoint-point).norm();
	}
	double c2 = v.squaredNorm();
	if (c2 <= c1)
	{
		closestPoint = segment.first;
		return (closestPoint - point).norm();
	}
	double b = c1 / c2;
	auto prb = segment.second + v*b;
	closestPoint = prb;
	return (prb-point).norm();
}

MyPolyline simplify(const MyPolyline & poly, double eps)
{
	if (poly.size() <= 3)
		return poly;

	std::vector<bool> mask(poly.size());
	std::queue<std::pair<int, int>> queue;
	// Find the point with the maximum distance
	queue.push({ 0,poly.size()-1});

	while (!queue.empty())
	{
		auto testPair = queue.front();
		queue.pop();

		double dmax = 0;
		int index = -1;

		for (int i = testPair.first+1; i < testPair.second; ++i)
		{
			double d = perpendicularDistance(poly[i], { poly[testPair.first],  poly[testPair.second] });
			if (d > dmax) {
				index = i;
				dmax = d;
			}
		}

		// If max distance is greater than epsilon, simplify
		if (dmax > eps) {
			queue.push({ testPair.first, index });
			queue.push({ index,testPair.second });
		}
		else {
			mask[testPair.first] = true;
			mask[testPair.second] = true;
		}

	}

	MyPolyline result;
	for (int i = 0; i < mask.size(); ++i)
	{
		if (mask[i])
			result.push_back(poly[i]);
	}
	return result;
}
