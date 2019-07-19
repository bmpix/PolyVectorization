#include "stdafx.h"
#include "typedefs.h"

Eigen::Vector2d _toEig(std::complex<double> c)
{
	 return Eigen::Vector2d(c.real(), c.imag());
}

Box findBoundingBox(const MyPolyline& poly)
{
	Box result;
	for (int i = 0; i < poly.size(); ++i)
	{
		result.xMin = std::min(result.xMin, poly[i].x());
		result.xMax = std::max(result.xMax, poly[i].x());
		result.yMin = std::min(result.yMin, poly[i].y());
		result.yMax = std::max(result.yMax, poly[i].y());
	}
	return result;
}

void bitwiseOr(BoolMatrix& lhs, const BoolMatrix& rhs)
{
#pragma omp parallel for
	for (int i = 0; i < lhs.rows(); ++i)
		for (int j = 0; j < lhs.cols(); ++j)
			lhs(i, j) = lhs(i, j) || rhs(i, j);
}

bool doesNeighborExist(int i, int j, int m, int n, bool vertical, bool left)
{
	int coord = vertical ? i : j;
	int bounds = vertical ? m : n;
	if (left)
		return coord > 0;
	else
		return coord < bounds - 1;
}

bool useNeighbor(int i, int j, int m, int n, bool vertical, bool left, const cv::Mat & mask, std::pair<int, int>& outNeighbor)
{
	int shift = left ? -1 : 1;
	int iNeigh = vertical ? i + shift : i;
	int jNeigh = !vertical ? j + shift : j;

	bool result = doesNeighborExist(i, j, m, n, vertical, left) && (mask.at<uchar>(iNeigh, jNeigh) != 0);
	if (result)
	{
		outNeighbor.first = iNeigh;
		outNeighbor.second = jNeigh;
		//outNeighbor = std::make_pair(iNeigh, jNeigh);
	}
	return result;
}

Eigen::Vector2d tangent(const MyPolyline& poly, int i)
{
	Eigen::Vector2d result = (i == poly.size() - 1 ? poly[i] - poly[i - 1] : poly[i + 1] - poly[i]);
	return result;
}

bool rootMatching(const std::array<std::complex<double>, 2>& r0, const std::array<std::complex<double>, 2>& r1)
{
	Eigen::Vector2d root00 = _toEig(r0[0]), root01 = _toEig(r0[1]), root10 = _toEig(r1[0]), root11 = _toEig(r1[1]);
	double dist0011 = std::max(fabs(root00.normalized().dot(root10.normalized())), fabs(root01.normalized().dot(root11.normalized())));
	double dist0110 = std::max(fabs(root00.normalized().dot(root11.normalized())), fabs(root01.normalized().dot(root10.normalized())));
	return dist0011 > dist0110;
}

bool rootMatching(const std::array<Eigen::MatrixXcd, 2>& roots, const std::array<int, 2>& p0, const std::array<int, 2>& p1)
{
	return rootMatching({ roots[0](p0[0], p0[1]), roots[1](p0[0], p0[1]) },
						{ roots[0](p1[0], p1[1]), roots[1](p1[0], p1[1]) });
}

bool inBounds(const Eigen::Vector2d & p, const cv::Mat & mask)
{
	int m = mask.rows, n = mask.cols;
	cv::Point cvPoint(std::round(p.x()), std::round(p.y()));
	return !((p.x() < 0) || (p.y() < 0) || (p.x() > n - 1) || (p.y() > m - 1) || (mask.at<uchar>(cvPoint) == 0));
}

bool inBounds(const std::array<int, 2>& p, const cv::Mat & mask)
{
	int m = mask.rows, n = mask.cols;
	return !((p[0] < 0) || (p[1] < 0) || (p[0] > m - 1) || (p[1] > n - 1) || (mask.at<uchar>(p[0],p[1]) == 0));
}

bool inBounds(int i, int j, const cv::Mat & mask)
{
	std::array<int, 2> p = { i,j };
	return inBounds(p, mask);
}
