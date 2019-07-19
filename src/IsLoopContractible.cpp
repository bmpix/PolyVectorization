#include "stdafx.h"
#include "IsLoopContractible.h"
#include "CircularSegment.h"
#include "Params.h"
std::vector<cv::Point2f> subcurve(const MyPolyline& poly, const CircularSegment& s)
{
	std::vector<cv::Point2f> real;
	Eigen::Vector2d start = poly[s.begin()];
	Eigen::Vector2d end = poly[s.end()];
	real.push_back({ (float)start.x(),(float)start.y() });
	for (int j = std::floor(s.begin()); j != std::floor(s.end()); j = (j + 1) % poly.size())
	{
		real.push_back({ (float)poly[j].x(),(float)poly[j].y() });
	}
	real.push_back({ (float)end.x(),(float)end.y() });
	return real;
}

const std::vector <cv::Point2f> convertToPolyline(const std::vector<edge_descriptor>& loop, const G& g, const std::vector<MyPolyline>& polys)
{
	std::vector<cv::Point2f> cvLoop;
	
	for (const auto& e : loop)
	{
		//now there's an edge bewteen two clusters
		//add the segment of the curve that connects them
		int myCurve = g[e].edgeCurve;

		/*std::cout << "Edge: " << e.m_source << " " << e.m_target << ", curve = " << myCurve << " " << "clPS: " << g[e.m_source].clusterPoints[0].curve << " @" << g[e.m_source].clusterPoints[0].segmentIdx << ", clPT: " <<
			g[e.m_target].clusterPoints[0].curve << " @" << g[e.m_target].clusterPoints[0].segmentIdx;
		if (myCurve != -1)
			std::cout << " (length = " << polys[myCurve].size() << ")";
		
		std::cout << std::endl;*/

		std::map<size_t, PointOnCurve> pt;
		std::vector<cv::Point2f> curve;
		if (myCurve != -1)
		{
			const auto& sourcePts = g[e.m_source].clusterPoints;
			const auto& targetPts = g[e.m_target].clusterPoints;
			CircularSegment bestS;
			double bestLength = std::numeric_limits<double>::max();
			bool curveClosed = (polys[myCurve].front() - polys[myCurve].back()).squaredNorm() < 1e-6;
			double modulus = curveClosed ? polys[myCurve].size() : -1.0;
			bool inverted;
			for (const auto& pt1 : sourcePts)
			{
				if (pt1.curve == myCurve)
				{
					for (const auto& pt2 : targetPts)
					{
						if (pt2.curve == myCurve)
						{
							CircularSegment s1 = CircularSegment(pt1.segmentIdx, pt2.segmentIdx, modulus);
							CircularSegment s2 = CircularSegment(pt2.segmentIdx, pt1.segmentIdx, modulus);
							double minLen = std::min(s1.length(), s2.length());
							if (minLen < bestLength)
							{
								if (s1.length() < s2.length())
								{
									bestS = s1;
									inverted = false;
								}
								else
								{
									bestS = s2;
									inverted = true;
								}
								bestLength = minLen;
							}
						}
					}
				}
			}

			curve = subcurve(polys[myCurve], bestS);

			if (inverted)
				std::reverse(curve.begin(), curve.end());

			/*
			int startPt = (pt[e.m_source].segmentIdx <= pt[e.m_target].segmentIdx) ? e.m_source : e.m_target;
			int endPt = (pt[e.m_source].segmentIdx > pt[e.m_target].segmentIdx) ? e.m_source : e.m_target;

			curve.push_back({ (float)pt[startPt].p.x(), (float)pt[startPt].p.y() });

			for (int i = std::floor(pt[startPt].segmentIdx + 1); i < pt[endPt].segmentIdx; ++i)
			{
			curve.push_back({ (float)polys[myCurve][i].x(), (float)polys[myCurve][i].y() });
			}

			cv::Point2f newPoint((float)pt[endPt].p.x(), (float)pt[endPt].p.y());
			if (!curve.empty() && (cv::norm(curve.back() - newPoint) > 1e-6))
			curve.push_back(newPoint);

			if (startPt != e.m_source)
			std::reverse(curve.begin(), curve.end());*/

		}
		else
		{
			curve = { { (float)g[e.m_source].location.x(),(float)g[e.m_source].location.y() },{ (float)g[e.m_target].location.x(),(float)g[e.m_target].location.y() } }; //for the cases around singularities
		}

		cvLoop.insert(cvLoop.end(), curve.begin(), curve.end());
	}
	return cvLoop;
}

bool isLoopContractible(const std::vector<edge_descriptor>& loop, const cv::Mat & origMask, const G& g, const std::vector<MyPolyline>& polys, std::vector<cv::Point2f>& cvLoop)
{
	//test if the loop contains any white pixels
	//essentially a scan conversion

	cvLoop = convertToPolyline(loop, g, polys);
	int minI = std::numeric_limits<int>::max(), maxI = std::numeric_limits<int>::min(),
		minJ = minI, maxJ = maxI;

	std::vector<cv::Point2f> cvLoopFiltered;
	for (int i = 0; i < cvLoop.size(); ++i)
	{
		if (cvLoopFiltered.empty() || (cv::norm(cvLoop[i] - cvLoopFiltered.back()) > 1e-6))
			cvLoopFiltered.push_back(cvLoop[i]);

		minI = std::min(minI, (int)cvLoop[i].y);
		minJ = std::min(minJ, (int)cvLoop[i].x);
		maxI = std::max(maxI, (int)cvLoop[i].y + 1);
		maxJ = std::max(maxJ, (int)cvLoop[i].x + 1);
	}
	cvLoop = cvLoopFiltered;
	if (!(minI >= 0 && minJ >= 0 && maxI <= origMask.rows && maxJ <= origMask.cols))
	{
		std::cout << "[isLoopContractible]: FATAL ERROR: " << minI << " " << minJ << " " << maxI << " " << maxJ << std::endl;
	}

	//if ((loop.size() > 10) && (loop.size() < 100))
	//{

	//	Box box;
	//	for (int i = 0; i < cvLoop.size(); ++i)
	//	{
	//		box.xMin = std::min(box.xMin, (double)cvLoop[i].x);
	//		box.xMax = std::max(box.xMax, (double)cvLoop[i].x);
	//		box.yMin = std::min(box.yMin, (double)cvLoop[i].y);
	//		box.yMax = std::max(box.yMax, (double)cvLoop[i].y);
	//	}
	//	
	//	svg::Document docRSPT("aloop.svg", svg::Layout(svg::Dimensions(box.xMax- box.xMin, box.yMax- box.yMin), svg::Layout::TopLeft));
	//	svg::Polyline poly(svg::Fill(), svg::Stroke(0.5, svg::Color::Black));
	//	for (size_t c = 0; c < cvLoop.size(); ++c)
	//	{
	//		poly << svg::Point(cvLoop[c].x, cvLoop[c].y);
	//	}
	//	docRSPT << poly;
	//	docRSPT.save();
	//}

	int numWhitePixels = 0;
	const int threshold = MAX_NUMBER_OF_WHITE_PIXELS_IN_A_CONTRACTIBLE_LOOP;
	for (int i = minI; i < maxI; ++i)
		for (int j = minJ; j < maxJ; ++j)
		{
			if (origMask.at<uchar>(i, j) == 0)
			{
				if (cv::pointPolygonTest(cvLoop, { (float)j,(float)i }, false) > 0)
				{
					++numWhitePixels;
					if (numWhitePixels > threshold)
						return false;
				}
			}
		}
	return (numWhitePixels <= threshold);
}
