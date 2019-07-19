#include "stdafx.h"
#include "AlmostReebGraph.h"
#include "intersections.h"
#include "simple_svg_1.0.0.hpp"
#include "chooseRoot.h"
#include "greedyTrace.h"
#include "ScanConvert.h"
#include "ChainDecomposition.h"
#include <ctime>

Cluster createCluster(int seedCurve, int seedPtIdx, const cv::Mat & origMask, const std::array<Eigen::MatrixXcd, 2>& roots, const std::vector<MyPolyline>& polys,
	const std::map<std::array<int, 2>, std::vector<PixelInfo>>& pixelInfo, const std::set<int>& onlyTheseCurves, const std::set<std::array<int, 2>>& singularities, const Eigen::MatrixXi& indices, const Eigen::VectorXcd& X)
{
	int m = origMask.rows, n = origMask.cols;
	Eigen::Vector2d myTangent = tangent(polys[seedCurve], seedPtIdx).normalized();
	Eigen::Vector2d pt = polys[seedCurve][seedPtIdx];
	int i = std::round(pt.y()), j = std::round(pt.x());

	int myRootIdx;
	double dot0 = fabs(myTangent.dot(_toEig(roots[0](i, j)).normalized())), dot1 = fabs(myTangent.dot(_toEig(roots[1](i, j)).normalized()));
	if (dot0 < dot1)
		myRootIdx = 1;
	else
		myRootIdx = 0;

	if ((fabs(dot0 - dot1) < 1e-1) || (singularities.find({ i,j }) != singularities.end()))
		return Cluster();

	Eigen::Vector2d traceDir(myTangent.y(), -myTangent.x());
	//shoot a perpendicular line
	std::map < std::array<int, 2>, bool > isMatchingOK;
	std::array<double, 2> closestSingularity = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
	//std::set<std::array<int, 2>> pixels = scanline(origMask, roots, pt, traceDir, isMatchingOK, bounds, singularities);
	std::map<std::array<int, 2>, std::vector<PixelInfo>> tmp1, tmp2;
	auto clusterCurve = greedyTrace(origMask, roots, pt, myTangent, X, tmp1, tmp2, indices, -1, true, closestSingularity).first;

	//find the initial point in the poly. totally unnecessary, but whatever
	int startIdx = -1;
	for (int i = 0; i < clusterCurve.size(); ++i)
	{
		if ((clusterCurve[i] - pt).squaredNorm() < 1e-10)
		{
			startIdx = i;
			break;
		}
	}
	assert(startIdx != -1);
	std::array<bool, 2> ignoreMe;
	std::set<std::array<int, 2>> pixels = scanConvert(clusterCurve, origMask, roots, startIdx, isMatchingOK, ignoreMe, singularities, { 0.0,(double)clusterCurve.size() - 1 });

	Eigen::ParametrizedLine<double, 2> pline(pt, traceDir);

	std::set<size_t> potentialCurves;
	for (const auto& p : pixels)
	{
		auto it = pixelInfo.find(p);
		if (it != pixelInfo.end())
		{
			for (int piIndex = 0; piIndex < it->second.size(); ++piIndex)
				potentialCurves.insert(it->second[piIndex].curve);
		}
	}

	Cluster result;

	//now make the intersections precise
	std::vector<std::pair<double, PointOnCurve>> curvesSortedByIntersectionPt;

	PointOnCurve seedPoint;
	seedPoint.curve = seedCurve;
	seedPoint.segmentIdx = seedPtIdx;
	seedPoint.p = pt;
	seedPoint.root = _toEig(roots[myRootIdx](i, j)).normalized();

	curvesSortedByIntersectionPt = { {0,seedPoint} };

	auto bbox = findBoundingBox(clusterCurve);

	for (size_t k : potentialCurves)
	{
		for (int segmentIdx = 0; segmentIdx + 1 < polys[k].size(); ++segmentIdx)
		{
			if ((k == seedCurve) && (fabs(segmentIdx - seedPtIdx) < 2))
				continue;

			if (bbox.completelyOutside(polys[k][segmentIdx],polys[k][segmentIdx+1]))
				continue;

			Box miniBox;
			std::tie(miniBox.xMin, miniBox.xMax) = std::minmax(polys[k][segmentIdx].x(), polys[k][segmentIdx + 1].x());
			std::tie(miniBox.yMin, miniBox.yMax) = std::minmax(polys[k][segmentIdx].y(), polys[k][segmentIdx + 1].y());

			for (int i = 0; i + 1 < clusterCurve.size(); ++i)
			{
				if (miniBox.completelyOutside(clusterCurve[i], clusterCurve[i + 1]))
					continue;

				double t, s;
				if (get_line_intersection(polys[k][segmentIdx].x(), polys[k][segmentIdx].y(), polys[k][segmentIdx + 1].x(), polys[k][segmentIdx + 1].y(),
					clusterCurve[i].x(), clusterCurve[i].y(), clusterCurve[i + 1].x(), clusterCurve[i + 1].y(), nullptr, nullptr, &s, &t))
				{
					PointOnCurve pc;
					Eigen::Vector2d edge = polys[k][segmentIdx + 1] - polys[k][segmentIdx];
					pc.p = polys[k][segmentIdx] + t*edge;
					pc.segmentIdx = segmentIdx + t;
					pc.curve = k;
					std::array<int, 2> pi = { std::round(pc.p.y()), std::round(pc.p.x()) };
					if (singularities.find(pi) != singularities.end()) //tough luck, skip it
						continue;

					if (origMask.at<uchar>(pi[0], pi[1]) == 0)
					{
						//there is that weird border case where the curve slightly goes outside the mask, and the intersection is right outside the mask
						bool fixed = false;
						const int width = 3;
						for (int i = -width / 2; i <= width / 2; ++i)
						{
							for (int j = -width / 2; j <= width / 2; ++j)
							{
								std::array<int, 2> newPixel = { pi[0] + i,pi[1] + j };
								if (!inBounds(newPixel, origMask) || (isMatchingOK.find(newPixel) == isMatchingOK.end()))
									continue;

								pi = newPixel;
								fixed = true;
								break;
							}
							if (fixed)
								break;
						}
						assert(fixed);
					}

					Eigen::Vector2d root0 = _toEig(roots[0](pi[0], pi[1])).normalized(), root1 = _toEig(roots[1](pi[0], pi[1])).normalized();
					int theirRootIdx = (fabs(edge.dot(root0)) < fabs(edge.dot(root1))) ? 1 : 0;
					pc.root = _toEig(roots[theirRootIdx](pi[0], pi[1])).normalized();

					assert(isMatchingOK.find(pi) != isMatchingOK.end());
					Eigen::Vector2d theirTangent = tangent(polys[pc.curve], pc.segmentIdx).normalized();
					if (onlyTheseCurves.empty() && ((isMatchingOK[pi] ^ (theirRootIdx == myRootIdx)) /*|| (fabs(theirTangent.normalized().dot(myTangent)) < 0.99)*/))
						continue;
					curvesSortedByIntersectionPt.push_back({ (i - startIdx + s)*0.1,pc });
				}

			}
		}
	}

	//now we have tons of intersection of AB with various curves
	//let's only keep my connected component of that
	std::sort(curvesSortedByIntersectionPt.begin(), curvesSortedByIntersectionPt.end(), [](const std::pair<double, PointOnCurve>& a, const std::pair<double, PointOnCurve>& b) {return a.first < b.first; });
	auto itK = std::find_if(curvesSortedByIntersectionPt.begin(), curvesSortedByIntersectionPt.end(), [&seedCurve](const std::pair<double, PointOnCurve>& a) {return (a.second.curve == seedCurve) && (fabs(a.first) < 1e-5); });
	assert(itK != curvesSortedByIntersectionPt.end());
	int kIdx = itK - curvesSortedByIntersectionPt.begin();

	std::array<bool, 2> hitASingularity = { fabs(closestSingularity[0] - curvesSortedByIntersectionPt[kIdx].first) < 1, fabs(closestSingularity[1] - curvesSortedByIntersectionPt[kIdx].first) < 1 };
	std::vector<PointOnCurve> myConnectedComponent = { seedPoint };
	int dirIdx = 0;
	for (int increment : {-1, 1})
	{
		for (int kk = kIdx + increment * 1; (kk < curvesSortedByIntersectionPt.size()) && (kk >= 0); kk += increment)
		{
			//if (fabs(curvesSortedByIntersectionPt[kk].first - curvesSortedByIntersectionPt[kk - increment].first) > 0.5*diam / (bounds[0] + bounds[1]))
			double diff = fabs(curvesSortedByIntersectionPt[kk].first - curvesSortedByIntersectionPt[kk - increment].first);
			if (fabs(closestSingularity[dirIdx] - curvesSortedByIntersectionPt[kk].first) < 1)
				hitASingularity[dirIdx] = true;

			if (diff > 1)
				break;
			else
			{
				if (increment == 1)
					myConnectedComponent.push_back(curvesSortedByIntersectionPt[kk].second);
				else
					myConnectedComponent.insert(myConnectedComponent.begin(), curvesSortedByIntersectionPt[kk].second);
			}
		}
		++dirIdx;
	}
	result.clusterPoints = myConnectedComponent;

	//select into a set all curves that this perpendicular intersects and that have the same root as me
	if (!result.clusterPoints.empty())
	{
		//orient all vectors
		Eigen::Vector2d vv = result.clusterPoints.begin()->root;
		for (auto& cp : result.clusterPoints)
		{
			if (vv.dot(cp.root) < 0)
				cp.root = -cp.root;
		}

		result.location = Eigen::Vector2d(0, 0);
		result.root = Eigen::Vector2d(0, 0);

		for (auto& cp : result.clusterPoints)
		{
			result.location += cp.p;
			result.root += cp.root;
		}

		result.location /= result.clusterPoints.size();
		//result.root /= result.clusterPoints.size();
		result.root = _toEig(roots[myRootIdx](i, j)).normalized();
		result.seedCurve = seedCurve;
		result.sharpCorner = false;
		result.nextToSingularity = false;
		result.clusterCurveHitSingularity = hitASingularity[0] | hitASingularity[1];
		result.width = (result.clusterPoints.back().p - result.clusterPoints.front().p).norm();
		result.clusterCurve = clusterCurve;
		result.split = false;
	}
	return result;
}

void contractSingularityBranches(G& g)
{
	std::map<edge_descriptor, size_t> myChain;
	auto chains = chainDecomposition(g, myChain);
	std::set<size_t> newSingularVertices;
	std::set<edge_descriptor> edgesToRemove;
	for (const auto& ch : chains)
	{
		if ((boost::degree(ch.front().m_source, g) == 1) ^ (boost::degree(ch.back().m_target, g) == 1))
		{
			auto chain = ch;
			if (boost::degree(ch.back().m_target, g) == 1)
				std::reverse(chain.begin(), chain.end()); //warning, directions of edges are screwed up

			for (const auto& e : chain)
			{
				if (g[e.m_source].clusterCurveHitSingularity || g[e.m_target].clusterCurveHitSingularity)
				{
					edgesToRemove.insert(e);
					newSingularVertices.insert(e.m_source);
					newSingularVertices.insert(e.m_target);
				}
				else
					break;
			}

		}
	}

	for (auto e : edgesToRemove)
		boost::remove_edge(e, g);

	for (size_t v : newSingularVertices)
		g[v].nextToSingularity = true;
}

void connectStuffAroundSingularities(G& g, const cv::Mat & origMask, const std::vector<MyPolyline>& polys, const std::set<std::array<int, 2>>& singularities, const std::array<Eigen::MatrixXcd, 2>& roots, const std::vector<std::array<bool, 2>>& endedWithASingularity)
{

	//add the extra singularities
	auto extendedSingularities = singularities;
	for (int i = 0; i < polys.size(); ++i)
	{
		if (endedWithASingularity[i][0])
			extendedSingularities.insert({ (int)polys[i].front().y(),(int)polys[i].front().x() });
		if (endedWithASingularity[i][1])
			extendedSingularities.insert({ (int)polys[i].back().y(),(int)polys[i].back().x() });
	}

	std::vector<std::pair<size_t, size_t>> newEdges;

	edge_iter eit, eend;
	for (std::tie(eit, eend) = boost::edges(g); eit != eend; ++eit)
		g[*eit].weight = 1.0;

	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
	

		if (g[v].nextToSingularity && (boost::degree(v, g) == 1))
		{
			//connect it to the closest singular vertex
			std::vector<vertex_descriptor> pDij(num_vertices(g));
			std::vector<double> dDij(num_vertices(g));
			auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx, g));
			auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx, g));
			dijkstra_shortest_paths(g, v,
				predecessor_map(predMap).
				distance_map(distMap).weight_map(get(&Edge::weight, g)));

			auto e = boost::out_edges(v, g).first;
			Eigen::Vector2d eVec = g[e->m_source].location - g[e->m_target].location;

			double bestDist = 10.0; //should be something reasonable
			int bestV = -1;
			for (size_t v1 = 0; v1 < boost::num_vertices(g); ++v1)
			{
				if (g[v1].nextToSingularity && !boost::edge(v, v1, g).second && (v1!=v) && boost::degree(v1,g)>=1)
				{
					//see if the straight line is within the narrow band and there is a singularity along the way
					std::vector<Eigen::Vector2d> poly;
					Eigen::Vector2d vec = g[v1].location - g[v].location;
					double dist = (g[v1].location - g[v].location).norm();
					//std::cout << "dist: " << dist << std::endl;
					int nSamples = std::max((int)(dist * 10),3);
					for (int j = 0; j < nSamples; ++j)
						poly.push_back(g[v].location+j*(vec / (nSamples - 1)));
					
					bool shouldAdd = true;
					for (int j = 0; j < poly.size(); ++j)
					{
						if (!inBounds(poly[j], origMask))
							shouldAdd = false;
					}

					if (shouldAdd && dDij[v1]>20 && vec.dot(eVec)>0)
					{
						std::map<std::array<int, 2>, bool> isRootMatchingOK;
						std::array<bool, 2> hitSingularity;
						scanConvert(poly, origMask, roots, 0, isRootMatchingOK, hitSingularity, extendedSingularities, { 0.0,(double)poly.size() - 1 });

						if (hitSingularity[0] | hitSingularity[1])
						{
							if (dist < bestDist)
							{
								bestV = v1;
								bestDist = dist;
							}
						}
					}
				}
			}


			if (bestV != -1)
			{
				newEdges.push_back({ v,(size_t)bestV });
				std::cout << "[connectStuffAroundSingularities]: Connecting vertex " << v << " to " << bestV << std::endl;
			}
		}
	}

	for (auto pair : newEdges)
	{
		auto e = boost::add_edge(pair.first, pair.second, g);
		g[e.first].edgeCurve = -1;
	}
}

G computeAlmostReebGraph(const cv::Mat & origMask, const std::array<Eigen::MatrixXcd, 2>& roots, const std::vector<MyPolyline>& polys, std::map<std::array<int, 2>, std::vector<PixelInfo>>& pixelInfo, const std::set<std::array<int, 2>>& singularities, const Eigen::MatrixXi& indices, const Eigen::VectorXcd& X, const std::vector<std::array<bool,2>>& endedWithASingularity)
{
	std::clock_t begin = -std::clock();
	std::cout << "Computing Reeb graph...";
	//int clusterIdx = 0;
	G resultList;
	typedef std::pair<G::vertex_descriptor, double> ClusterAndItsIntersection;
	std::map<int, std::vector<ClusterAndItsIntersection>> intersectionsPerCurve;


	auto areClustersEquivalent = [](const Cluster& a, const Cluster& b)
	{
		if (a.clusterPoints.size() != b.clusterPoints.size())
			return false;
		for (int i = 0; i < a.clusterPoints.size(); ++i)
			if (a.clusterPoints[i].curve != b.clusterPoints[i].curve)
				return false;
		return true;
	};

	std::map < std::array<int, 2>, std::set<int>> clustersPerPixel;
	std::cout << std::endl;

	std::vector<PointOnCurve> seedPts;
	for (int k = 0; k < polys.size(); ++k)
	{
		PointOnCurve pc;
		pc.curve = k;
		pc.segmentIdx = 0;
		seedPts.push_back(pc);
		pc.segmentIdx = polys[k].size() - 1;
		seedPts.push_back(pc);
	}

	std::vector<Cluster> clusters(seedPts.size());
	//for (int k = 0; k < polys.size(); ++k)
#pragma omp parallel for
	for (int clusterIdx = 0; clusterIdx < seedPts.size(); ++clusterIdx)
	{
		Cluster cl = createCluster(seedPts[clusterIdx].curve, seedPts[clusterIdx].segmentIdx, origMask, roots, polys, pixelInfo, {}, singularities, indices, X);

		std::array<int, 2> p0i = { (int)std::round(cl.location.y()),(int)std::round(cl.location.x()) };
		cl.clusterIdx = clusterIdx;
		clusters[clusterIdx] = cl;
	}

	for (int i = 0; i < clusters.size(); ++i)
	{
		auto v = boost::add_vertex(resultList);
		resultList[v] = clusters[i];

		for (int j = 0; j < clusters[i].clusterPoints.size(); ++j)
			intersectionsPerCurve[clusters[i].clusterPoints[j].curve].push_back(ClusterAndItsIntersection(v, clusters[i].clusterPoints[j].segmentIdx));
	}
	//now all clusters are ready
	//let's order them
	for (int i = 0; i < polys.size(); ++i)
	{
		std::sort(intersectionsPerCurve[i].begin(), intersectionsPerCurve[i].end(), [](const ClusterAndItsIntersection& a, const ClusterAndItsIntersection& b) {return a.second < b.second; });
		for (int j = 0; j < static_cast<int>(intersectionsPerCurve[i].size()) - 1; ++j)
		{
			auto pair = std::make_pair(intersectionsPerCurve[i][j].first, intersectionsPerCurve[i][j + 1].first);
			if (!boost::edge(pair.first, pair.second, resultList).second) //avoid duplicate edges
			{
				auto edge = boost::add_edge(pair.first, pair.second, resultList);
				resultList[edge.first].edgeCurve = i;
			}
		}
	}

	for (int i = 0; i < polys.size(); ++i)
	{
		if (!intersectionsPerCurve[i].empty() && (endedWithASingularity[i][0] | endedWithASingularity[i][1]))
		{
			if (endedWithASingularity[i][0])
				resultList[intersectionsPerCurve[i].front().first].nextToSingularity = true;
			if (endedWithASingularity[i][1])
				resultList[intersectionsPerCurve[i].back().first].nextToSingularity  = true;
		}
	}

	std::cout << "done in " << double(begin + clock()) / CLOCKS_PER_SEC << " seconds. " << std::endl;
	return resultList;
}

void simpleThresholds(G& g)
{
	//simple thresholding: if it's a single curve => delete
	std::vector<int> components(boost::num_vertices(g));
	int nComponents = boost::connected_components(g, &components[0]);
	std::vector<Box> bboxes(nComponents);
	std::vector<bool> shouldntRemove(nComponents);
	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		if (g[v].clusterPoints.size() > 1)
			shouldntRemove[components[v]] = true;

		bboxes[components[v]].xMax = std::max(bboxes[components[v]].xMax, g[v].location.x());
		bboxes[components[v]].xMin = std::min(bboxes[components[v]].xMin, g[v].location.x());
		bboxes[components[v]].yMax = std::max(bboxes[components[v]].yMax, g[v].location.y());
		bboxes[components[v]].yMin = std::min(bboxes[components[v]].yMin, g[v].location.y());
	}
	
	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		const auto& bbox = bboxes[components[v]];
		if ((!shouldntRemove[components[v]]) || ((bbox.xMax-bbox.xMin)*(bbox.yMax-bbox.yMin)<10))
			boost::clear_vertex(v, g);
	}
}
