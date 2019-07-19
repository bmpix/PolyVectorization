#include "chopFakeEnds.h"
#include "greedyTrace.h"
#include "intersections.h"
#include <iostream>
#include "Params.h"
#include "graph_typedefs.h"
#include "ContractLoops.h"
std::pair<std::vector<MyPolyline>, G> chopFakeEnds(const std::vector<MyPolyline>& polys, const std::vector<std::vector<double>>& radii, const std::vector<std::array<bool, 2>>& protectedEnds,
	const std::vector<std::array<bool, 2>>& isItASpecialDeg2Vertex, const std::vector<std::pair<PointOnCurve, PointOnCurve>>& yJunctions)
{
	std::vector<MyPolyline> result(polys.size());

	auto covered = [&polys, &radii](int i, double startSegment, double endSegment, int j)
	{
		double start = std::min(startSegment, endSegment);
		double end = std::max(startSegment, endSegment);

		double coveredLength = 0, totalLength = 0;
		for (int k = start; k < end; ++k)
		{
			bool vertexCovered = false;
			for (int k1 = 0; k1 < polys[j].size(); ++k1)
			{
				if ((polys[j][k1] - polys[i][k]).squaredNorm() < std::pow(radii[i][k] + radii[j][k1], 2))
				{
					vertexCovered = true;
					break;
				}
			}

			double length = (polys[i][k] - polys[i][k + 1]).norm();
			if (vertexCovered)
				coveredLength += length;

			totalLength += length;
		}

		if (((totalLength - coveredLength > 1) && (coveredLength / totalLength < PRUNE_SHORT_BRANCHES_RATIO)) || (totalLength > 10))
			return false;
		return true;
	};

	std::vector<Box> bboxes;
	for (int i = 0; i < polys.size(); ++i)
		bboxes.push_back(findBoundingBox(polys[i]));

	std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> extraEdges;
	std::vector<std::array<std::pair<double, PointOnCurve>, 2>> limits(polys.size());

	std::vector<std::array<bool, 2>> isItASpecialDef2VertexUpdated = isItASpecialDeg2Vertex; //if we cut this endpoint, it's no longer special

	for (int i = 0; i < polys.size(); ++i)
	{
		limits[i][0].first = 0;
		limits[i][0].second.curve = -1;
		limits[i][1].first = polys[i].size();
		limits[i][1].second.curve = -1;

		if (polys[i].size() <= 3) //ignore the short guys
			continue;

		for (int endIdx = 0; endIdx < 2; ++endIdx)
		{
			int startSegment = endIdx == 0 ? 0 : polys[i].size() - 1;
			int endSegment = polys[i].size() / 2;
			int incSegment = endIdx == 0 ? 1 : -1;

			bool foundLimit = false;
			for (int segmentIdx = startSegment; segmentIdx != endSegment; segmentIdx += incSegment)
			{
				Eigen::Vector2d s0 = polys[i][segmentIdx], s1 = polys[i][segmentIdx + incSegment];
				//now let's intersect
				typedef std::pair<double, PointOnCurve> A;
				std::vector<A> intersections;

				for (int j = 0; j < polys.size(); ++j)
				{
					//if (i == j) //todo: remove this
					//	continue;

					if (bboxes[j].completelyOutside(s0, s1))
						continue;

					for (int jSegmentIdx = 0; jSegmentIdx + 1 < polys[j].size(); ++jSegmentIdx)
					{
						if ((i == j) && (fabs(jSegmentIdx - segmentIdx) < 20)) //skip imaginary self-intersections, 20 is a random number that is big enough. can be resolved by careful engineering
							continue;

						Eigen::Vector2d sj0 = polys[j][jSegmentIdx], sj1 = polys[j][jSegmentIdx + 1];
						if (bboxes[i].completelyOutside(sj0, sj1))
							continue;

						double s, t;
						if (get_line_intersection(s0.x(), s0.y(), s1.x(), s1.y(), sj0.x(), sj0.y(), sj1.x(), sj1.y(), nullptr, nullptr, &s, &t))
						{
							if (endIdx == 1)
								t = 1 - t;

							if ((jSegmentIdx != 0 || s > 1e-6) && (jSegmentIdx != polys[j].size() - 1 || s < 1 - 1e-6) && (t > 1e-6) && (t < 1 - 1e-6) && (std::min(segmentIdx, segmentIdx + incSegment) + t > 1e-6))
							{
								PointOnCurve pt;
								pt.curve = j;
								pt.segmentIdx = jSegmentIdx + s;
								intersections.push_back({ std::min(segmentIdx,segmentIdx + incSegment) + t,pt });
							}
						}
					}
				}

				std::sort(intersections.begin(), intersections.end(), [](const A& a, const A& b) {return a.first < b.first; });
				if (endIdx == 1)
					std::reverse(intersections.begin(), intersections.end());

				if (!intersections.empty())
				{
					std::cout << "Curve " << i << " intersections (endIdx = " << endIdx << "): " << std::endl;
					for (auto ii : intersections)
						std::cout << "c" << ii.second.curve << " at " << ii.first << "(@" << ii.second.segmentIdx << ")" << std::endl;
				}

				for (int k = 0; k < intersections.size(); ++k)
				{
					if (!protectedEnds[i][endIdx] && !foundLimit && covered(i, startSegment, intersections[k].first, intersections[k].second.curve))
						limits[i][endIdx] = intersections[k];

					foundLimit = true;

					if (i != intersections[k].second.curve) //remove self-loops
						extraEdges.push_back(std::make_pair(std::make_pair(i, intersections[k].first), std::make_pair(intersections[k].second.curve, std::floor(intersections[k].second.segmentIdx))));
				}
			}
		}

		std::cout << "Curve " << i << " limits: " << limits[i][0].first << " and " << limits[i][1].first << "(length: " << polys[i].size() << ")" << std::endl;
		MyPolyline newPoly;
		int start = limits[i][0].first > 1e-5 ? std::ceil(limits[i][0].first) : 0;
		int end = limits[i][1].first < polys[i].size() - 1e-5 ? std::floor(limits[i][1].first) + 1 : polys[i].size();
		newPoly.insert(newPoly.end(), polys[i].begin() + start, polys[i].begin() + end);

		if (newPoly.size() >= 2)
		{
			if (limits[i][0].first > 1e-5)
			{
				double lAlpha = limits[i][0].first - std::floor(limits[i][0].first);
				Eigen::Vector2d pCeil = polys[i][std::ceil(limits[i][0].first)];
				Eigen::Vector2d pFloor = polys[i][std::floor(limits[i][0].first)];
				newPoly.front() = pFloor + lAlpha*(pCeil - pFloor);
				isItASpecialDef2VertexUpdated[i][0] = false;
			}
			if (limits[i][1].first < polys[i].size() - 1e-5)
			{
				double lAlpha = limits[i][1].first - std::floor(limits[i][1].first);
				Eigen::Vector2d pCeil = polys[i][std::ceil(limits[i][1].first)];
				Eigen::Vector2d pFloor = polys[i][std::floor(limits[i][1].first)];
				newPoly.back() = pFloor + lAlpha*(pCeil - pFloor);
				isItASpecialDef2VertexUpdated[i][1] = false;
			}
		}
		result[i] = newPoly;
	}
	std::cout << "Before adding edges from Y-junctions: " << extraEdges.size() << std::endl;
	//add extraEdges from yJunctions
	for (int i = 0; i < yJunctions.size(); ++i)
	{
		int curve1 = yJunctions[i].first.curve, curve2 = yJunctions[i].second.curve;
		if ((limits[curve1][1].first <= yJunctions[i].first.segmentIdx) || (yJunctions[i].first.segmentIdx < limits[curve1][0].first))
			continue;

		if ((limits[curve2][1].first <= yJunctions[i].second.segmentIdx) || (yJunctions[i].second.segmentIdx < limits[curve2][0].first))
			continue;

		extraEdges.push_back(std::make_pair(std::make_pair(yJunctions[i].first.curve, yJunctions[i].first.segmentIdx), std::make_pair(yJunctions[i].second.curve, yJunctions[i].second.segmentIdx)));
	}

	//adjust the extra edges segment indices because of the chopped beginning of the curve
	for (int i = 0; i < extraEdges.size(); ++i)
	{
		int curve1 = extraEdges[i].first.first;
		//std::cout << "Curve " << curve1 << " adjusted int pt from " << extraEdges[i].first.second << " to ";
		extraEdges[i].first.second -= limits[curve1][0].first;
		extraEdges[i].first.second = std::max(extraEdges[i].first.second, 0);
		//std::cout << extraEdges[i].first.second << std::endl;

		int curve2 = extraEdges[i].second.first;
		//std::cout << "Curve " << curve2 << " adjusted int pt from " << extraEdges[i].second.second << " to ";
		extraEdges[i].second.second -= limits[curve2][0].first;
		extraEdges[i].second.second = std::max(extraEdges[i].second.second, 0);
		//std::cout << extraEdges[i].second.second << std::endl;
	}

	/*std::cout << "Extra edges: " << std::endl;
	for (int i = 0; i < extraEdges.size(); ++i)
	std::cout << extraEdges[i].first.first << "  " << extraEdges[i].first.second << " " << extraEdges[i].second.first << " " << extraEdges[i].second.second << std::endl;*/

	std::map<std::pair<int, int>, int> curveAndPtToIdx;

	G curveGraph;
	std::vector<std::vector<std::pair<double, int>>> curveVertices(result.size()); //for each curve, vertices in the new graph will be endpoints and all the intersections needed
	for (int i = 0; i < result.size(); ++i)
	{
		if (result[i].empty())
			continue;

		size_t v1 = boost::add_vertex(curveGraph);
		PointOnCurve pt1, pt2;
		pt1.curve = i; pt2.curve = i;
		pt1.segmentIdx = 0; pt2.segmentIdx = result[i].size() - 1;
		curveGraph[v1].clusterPoints = { pt1 };
		curveGraph[v1].location = result[i].front();
		curveVertices[i].push_back(std::make_pair(0.0, v1));
		curveAndPtToIdx[std::make_pair(i, 0)] = v1;
		curveGraph[v1].split = isItASpecialDef2VertexUpdated[i][0];

		size_t v2 = boost::add_vertex(curveGraph);
		curveGraph[v2].clusterPoints = { pt2 };
		curveGraph[v2].location = result[i].back();
		curveVertices[i].push_back(std::make_pair(result[i].size() - 1, v2));
		curveAndPtToIdx[std::make_pair(i, (int)(result[i].size() - 1))] = v2;
		curveGraph[v2].split = isItASpecialDef2VertexUpdated[i][1];
	}

	std::cout << "Orig vertices: " << boost::num_vertices(curveGraph) << std::endl;

	for (int jj = 0; jj < extraEdges.size(); ++jj)
	{
		for (const auto& intersectionPt : { extraEdges[jj].first, extraEdges[jj].second })
		{
			int curve = intersectionPt.first;
			int roundedOff = std::floor(intersectionPt.second);
			if (intersectionPt.second == 0 || intersectionPt.second == result[curve].size() - 1 || curveAndPtToIdx.find(std::make_pair(curve, roundedOff)) != curveAndPtToIdx.end())
				continue; //nothing to add

			if (result[curve].size() <= intersectionPt.second)
			{
				std::cout << "Skipping an extra edge connecting to curve " << curve << " at " << intersectionPt.second << "(length = " << result[curve].size() << ") " << std::endl;
				continue;
			}

			size_t v = boost::add_vertex(curveGraph);
			curveGraph[v].split = false;
			PointOnCurve pt;
			pt.curve = curve;
			pt.segmentIdx = intersectionPt.second;
			pt.p = result[curve][roundedOff];

			curveGraph[v].clusterPoints = { pt };
			curveGraph[v].location = result[curve][roundedOff];
			curveVertices[curve].push_back(std::make_pair(pt.segmentIdx, v));
			curveAndPtToIdx[std::make_pair(curve, roundedOff)] = v;
		}
	}

	for (int i = 0; i < boost::num_vertices(curveGraph); ++i)
	{
		curveGraph[i].clusterIdx = i;
		curveGraph[i].clusterCurveHitSingularity = false;
		curveGraph[i].nextToSingularity = false;
		curveGraph[i].root = Eigen::Vector2d(0, 1);
		curveGraph[i].seedCurve = 0;
		curveGraph[i].sharpCorner = false;
		curveGraph[i].width = 0;
	}

	for (int i = 0; i < boost::num_vertices(curveGraph); ++i)
	{
		std::cout << "Vertex " << i << ": " << curveGraph[i].clusterPoints[0].curve << " @" << curveGraph[i].clusterPoints[0].segmentIdx << "(length = " << result[curveGraph[i].clusterPoints[0].curve].size() << ")" << std::endl;
	}

	for (int i = 0; i < result.size(); ++i)
	{
		std::cout << "Curve " << i << " ";
		std::sort(curveVertices[i].begin(), curveVertices[i].end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {return a.first < b.first; });
		for (auto& jj : curveVertices[i])
		{
			std::cout << jj.second << " (@" << jj.first << ") ";
		}
		std::cout << std::endl;

		for (int j = 0; j + 1 < curveVertices[i].size(); ++j)
		{
			if (curveVertices[i][j].second == curveVertices[i][j + 1].second)
				continue;

			auto e = boost::add_edge(curveVertices[i][j].second, curveVertices[i][j + 1].second, curveGraph);
			curveGraph[e.first].edgeCurve = i;
			std::cout << "E " << curveVertices[i][j].second << " " << curveVertices[i][j + 1].second << std::endl;
		}
	}

	std::map<std::pair<size_t, size_t>, bool> added;

	for (int i = 0; i < extraEdges.size(); ++i)
	{
		extraEdges[i].first.second = std::floor(extraEdges[i].first.second);
		extraEdges[i].second.second = std::floor(extraEdges[i].second.second);
		if (curveAndPtToIdx.find(extraEdges[i].first) != curveAndPtToIdx.end() && curveAndPtToIdx.find(extraEdges[i].second) != curveAndPtToIdx.end())
		{
			if (!added[std::minmax(curveAndPtToIdx[extraEdges[i].first], curveAndPtToIdx[extraEdges[i].second])])
			{
				added[std::minmax(curveAndPtToIdx[extraEdges[i].first], curveAndPtToIdx[extraEdges[i].second])] = true;
				auto e = boost::add_edge(curveAndPtToIdx[extraEdges[i].first], curveAndPtToIdx[extraEdges[i].second], curveGraph);
				curveGraph[e.first].edgeCurve = -1;
				std::cout << "Adding extra edge: " << curveAndPtToIdx[extraEdges[i].first] << " " << curveAndPtToIdx[extraEdges[i].second] << std::endl;
			}
		}
	}

	return std::make_pair(result, curveGraph);

}
