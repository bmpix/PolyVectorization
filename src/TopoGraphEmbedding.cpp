#include "stdafx.h"
#include "TopoGraphEmbedding.h"
#include "ContractDeg2.h"
#include <ctime>

std::map<size_t, size_t> Distances::orient(const G& g, const std::vector<int>& chains, const std::set<size_t>& vertSet, int root, const std::set<edge_descriptor>& cutThoseEdges)
{
	//choose a root
	MyWeirdGraph subgraph;
	std::map<size_t, size_t> globalToLocal, localToGlobal;
	for (auto v : vertSet)
	{
		size_t i = boost::add_vertex(subgraph);
		globalToLocal[v] = i;
		localToGlobal[i] = v;
	}
	auto wMap = get(boost::edge_weight, subgraph);
	for (int chainIdx : chains)
	{
		for (const auto& e : chainsSeparated[chainIdx])
		{
			if (cutThoseEdges.find(e) == cutThoseEdges.end())
			{
				auto e1 = boost::add_edge(globalToLocal[e.m_source], globalToLocal[e.m_target], subgraph);
				auto e2 = boost::add_edge(globalToLocal[e.m_target], globalToLocal[e.m_source], subgraph);
				wMap[e1.first] = 1;
				wMap[e2.first] = 1;
			}
		}
	}

	std::vector<vertex_descriptor> parents(vertSet.size());
	std::vector<double> distances(vertSet.size());
	auto indexMap = boost::get(boost::vertex_index, subgraph);
	auto predMap = boost::make_iterator_property_map(parents.begin(), indexMap);
	auto distMap = boost::make_iterator_property_map(distances.begin(), indexMap);

	boost::dijkstra_shortest_paths(subgraph, globalToLocal[root], predecessor_map(predMap).distance_map(distMap));
	std::map<size_t, size_t> result;
	for (int i = 0; i < parents.size(); ++i)
		result[localToGlobal[i]] = localToGlobal[parents[i]];
	return result;
}

std::pair<double, std::pair<PointOnCurve, PointOnCurve>> Distances::distanceBetweenSamples(int v1, int s1, int v2, int s2)
{
	PointOnCurve bestPt1, bestPt2;
	double bestDist = 1e11;
	for (int curve : pairOfVerticesToSharedCurves[std::minmax(v1, v2)])
	{
		const auto& is = vertexAndCurveToSamples[std::make_pair(v1, curve)];
		const auto& js = vertexAndCurveToSamples[std::make_pair(v2, curve)];
		for (int i : is)
		{
			double iDist = (g[v1].clusterPoints[i].p - g[v1].clusterPoints[s1].p).norm();
			for (int j : js)
			{
				double dist = iDist + (g[v2].clusterPoints[j].p - g[v2].clusterPoints[s2].p).norm();

				if (fabs(g[v1].clusterPoints[i].segmentIdx - g[v2].clusterPoints[j].segmentIdx) > 500) //this is a hack to fight off spiraling curves
					dist = std::numeric_limits<double>::max();

				if (bestDist > dist)
				{
					bestDist = dist;
					bestPt1 = g[v1].clusterPoints[i];
					bestPt2 = g[v2].clusterPoints[j];
				}
			}
		}
	}

	if ((bestDist > 1e10) && (bestDist < 1e20)) //no shared curve found. this can only happen for edges connecting 'singular' vertices
	{
		double bestV1 = bestDist, bestV2 = bestDist;
		for (int i = 0; i < g[v1].clusterPoints.size(); ++i)
		{
			double dist = (g[v1].clusterPoints[i].p - g[v1].location).norm();
			if (dist < bestV1)
			{
				bestPt1 = g[v1].clusterPoints[i];
				bestV1 = dist;
			}
		}
		for (int i = 0; i < g[v2].clusterPoints.size(); ++i)
		{
			double dist = (g[v2].clusterPoints[i].p - g[v2].location).norm();
			if (dist < bestV2)
			{
				bestPt2 = g[v2].clusterPoints[i];
				bestV2 = dist;
			}
		}
		bestDist = (bestPt1.p - g[v1].clusterPoints[s1].p).norm() + (bestPt2.p - g[v2].clusterPoints[s2].p).norm();
	}

	return{ bestDist,{ bestPt1, bestPt2 } };
}

void Distances::fillInSeedTypes()
{
	int nSeeds = seeds.size();
	seedType.resize(nSeeds, ST_REGULAR);
	for (int seedIdx = 0; seedIdx < nSeeds; ++seedIdx)
	{
		auto seed = seeds[seedIdx];
		if (boost::degree(seed, g) <= 2)
			seedType[seedIdx] = ST_VAL2;
		else
		{
			//test if the seed has a loop adjacent
			std::map<int, int> valence;
			for (auto j : adjChains[seed])
			{
				valence[chainsSeparated[j].front().m_source]++;
				valence[chainsSeparated[j].back().m_target]++;
			}
			valence[seed] = 0;
			for (auto it : valence)
			{
				if ((it.second > 0) && (boost::degree(it.first, g) > 2))
					seedType[seedIdx] = ST_ADJ_TO_ANOTHER_SEED;
				else if (it.second > 1)
					seedType[seedIdx] = ST_LOOP;
			}
		}
	}
}

Distances::Distances(const G & g, const std::vector<std::pair<size_t, size_t>>& topoGraph, const std::vector<std::vector<edge_descriptor>>& chainsSeparated, const cv::Mat& bwImg)
	:g(g), topoGraph(topoGraph), chainsSeparated(chainsSeparated)
{
	std::clock_t begin = -std::clock();
	std::cout << "Computing lots of distances.." << std::endl;

	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		for (int i = 0; i < g[v].clusterPoints.size(); ++i)
		{
			vertexAndCurveToSamples[std::make_pair(v, g[v].clusterPoints[i].curve)].push_back(i);
		}
	}

	edge_iter eit, eend;
	for (std::tie(eit, eend) = boost::edges(g); eit != eend; ++eit)
	{
		std::set<int> curvesSource, curvesTarget;
		for (const auto& v : g[eit->m_source].clusterPoints)
			curvesSource.insert(v.curve);
		for (const auto& v : g[eit->m_target].clusterPoints)
			curvesTarget.insert(v.curve);

		std::vector<int> setIntersection(std::max(curvesSource.size(), curvesTarget.size()));
		auto it = std::set_intersection(curvesSource.begin(), curvesSource.end(), curvesTarget.begin(), curvesTarget.end(), setIntersection.begin());
		setIntersection.resize(it - setIntersection.begin());
		pairOfVerticesToSharedCurves[std::minmax(eit->m_source, eit->m_target)] = setIntersection;
	}

	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		if (boost::degree(v, g) > 2)
			seeds.push_back(v);
	}

	std::vector<bool> chainTaken(chainsSeparated.size());
	for (size_t seed : seeds)
	{
		//collect all edges and vertices of the component
		for (int j = 0; j < topoGraph.size(); ++j)
		{
			if (topoGraph[j].first == seed ^ topoGraph[j].second == seed)
			{
				adjChains[seed].push_back(j);
				chainTaken[j] = true;
			}
		}
	}

	for (int j = 0; j < chainsSeparated.size(); ++j)
		if (!chainTaken[j])
		{
			size_t seed = chainsSeparated[j].front().m_source;
			seeds.push_back(seed);
			adjChains[seed] = { j };
		}

	int nSeeds = seeds.size();

	graphs.resize(nSeeds);
	indexToClusterAndPt.resize(nSeeds);
	clusterAndPtToIndex.resize(nSeeds);
	vertexSets.resize(nSeeds);

	D.resize(nSeeds);
	parent.resize(nSeeds);

	fillInSeedTypes();

	std::vector<int> nSamples(nSeeds);
	std::vector<std::vector<std::pair<size_t, size_t>>> allEdges(nSeeds);
	std::vector<std::vector<double>> allWeights(nSeeds);
#pragma omp parallel for
	for (int seedIdx = 0; seedIdx < nSeeds; ++seedIdx)
	{
		size_t seed = seeds[seedIdx];

		//ENUMERATE ALL VERTICES ADJACENT TO THE SEED
		std::set<size_t> vertSet;
		for (const auto& ch : adjChains[seed])
			for (const auto& e : chainsSeparated[ch])
			{
				for (size_t v : {e.m_source, e.m_target})
					vertSet.insert(v);
			}

		vertexSets[seedIdx] = vertSet;

		//enumerate all possible locations
		int N = 0;
		for (size_t v : vertSet)
		{
			for (int i = 0; i < g[v].clusterPoints.size(); ++i)
			{
				clusterAndPtToIndex[seedIdx][std::make_pair(v, i)] = N;
				indexToClusterAndPt[seedIdx][N] = std::make_pair(v, i);
				++N;

			}
		}
		nSamples[seedIdx] = N;
		//std::cout << "Creating graph with " << N << " vertices ";

		//gather all edges of the distance graph
		auto& edges = allEdges[seedIdx];
		auto& weights = allWeights[seedIdx];

		for (const auto& ch : adjChains[seed])
		{
			std::vector<double> partialLengths(chainsSeparated[ch].size());
			partialLengths[0] = fabs((g[chainsSeparated[ch][0].m_target].location - g[chainsSeparated[ch][0].m_source].location).dot(g[chainsSeparated[ch][0].m_source].root));
			for (int k = 1; k < chainsSeparated[ch].size(); ++k)
			{
				const auto& e = chainsSeparated[ch][k];
				partialLengths[k] = partialLengths[k - 1] + fabs((g[e.m_target].location - g[e.m_source].location).dot(g[e.m_source].root));
			}

			for (int k = 0; k < chainsSeparated[ch].size(); ++k)
			{
				const auto& e = chainsSeparated[ch][k];
				double extraW = 1.0;
				if ((boost::degree(chainsSeparated[ch].front().m_source, g) == 1) && (partialLengths[k] < 1))
					extraW = 1e-9;
				if ((boost::degree(chainsSeparated[ch].back().m_target, g) == 1) && (k > 0) && (partialLengths.back() - partialLengths[k - 1] < 1))
					extraW = 1e-9;

				double strokeWidth1 = std::max(1.0, (g[e.m_source].clusterPoints.back().p - g[e.m_source].clusterPoints.front().p).norm());
				double strokeWidth2 = std::max(1.0, (g[e.m_target].clusterPoints.back().p - g[e.m_target].clusterPoints.front().p).norm());

				for (int i = 0; i < g[e.m_source].clusterPoints.size(); ++i)
				{
					double centerShiftI = (g[e.m_source].location - g[e.m_source].clusterPoints[i].p).norm();

					for (int j = 0; j < g[e.m_target].clusterPoints.size(); ++j)
					{
						double centerShiftJ = (g[e.m_target].location - g[e.m_target].clusterPoints[j].p).norm();

						edges.push_back(std::make_pair(clusterAndPtToIndex[seedIdx][std::make_pair(e.m_source, i)], clusterAndPtToIndex[seedIdx][std::make_pair(e.m_target, j)]));
						auto distAndInfo = distanceBetweenSamples(e.m_source, i, e.m_target, j);

						double centeringWeight = 0.07*(centerShiftI + centerShiftJ) / (strokeWidth1 + strokeWidth2);

						if (g[e.m_target].split && g[e.m_source].split)
							centeringWeight = 0;

						weights.push_back(extraW*(distAndInfo.first + centeringWeight));
					}
				}
			}
		}
	}

	//shouldn't be parallelized
	for (int seedIdx = 0; seedIdx < nSeeds; ++seedIdx)
	{
		size_t seed = seeds[seedIdx];
		const auto& edges = allEdges[seedIdx];
		const auto& weights = allWeights[seedIdx];
		if (edges.empty())
			continue;

		//std::cout << " and " << edges.size() << " edges" << std::endl;
		int N = nSamples[seedIdx];

		std::vector <int> d(N, (std::numeric_limits <double>::max)());

		for (int j : adjChains[seed])
		{
			size_t v;
			//distances should be computed from the fixed vertices
			if (seedType[seedIdx] != ST_VAL2)
				v = topoGraph[j].first == seed ? topoGraph[j].second : topoGraph[j].first;
			else
				v = seed;

			//orient all edges starting from that vertex

			//take an edge on the other side of the shared vertex
			//mark its weights as infinity
			oedge_iter eit, eend;
			std::set<MyWeirdGraph::edge_descriptor> cutThoseEdges_sampleGraph;
			std::set<edge_descriptor> cutThoseEdges_vertexGraph;
			boost::graph_traits < MyWeirdGraph >::edge_iterator e, e_end;
			for (std::tie(eit, eend) = boost::out_edges(v, g); eit != eend; ++eit)
			{
				if (std::find(chainsSeparated[j].begin(), chainsSeparated[j].end(), *eit) == chainsSeparated[j].end()) //there's an edge adjacent to v that's not in my chain
				{
					cutThoseEdges_vertexGraph.insert(*eit);
				}
			}

			auto vertexParent = orient(g, adjChains[seed], vertexSets[seedIdx], v, cutThoseEdges_vertexGraph);
			auto orientedEdges = edges;
			for (auto& e : orientedEdges)
			{
				int vtx1 = indexToClusterAndPt[seedIdx][e.first].first, vtx2 = indexToClusterAndPt[seedIdx][e.second].first;
				if ((vtx1 != vtx2) && (vertexParent[vtx2] != vtx1))
					e = { e.second, e.first };
			}

			printf("Creating graph with %d edges\r", orientedEdges.size());

			graphs[seedIdx][v] = std::make_shared<MyWeirdGraph>(&orientedEdges[0], orientedEdges.size() + &orientedEdges[0], N);

			for (const auto& eit : cutThoseEdges_vertexGraph)
			{
				for (std::tie(e, e_end) = boost::edges(*graphs[seedIdx][v]); e != e_end; ++e)
				{
					int index1 = indexToClusterAndPt[seedIdx][e->m_source].first, index2 = indexToClusterAndPt[seedIdx][e->m_target].first;
					if ((index1*index2 == eit.m_source*eit.m_target) && (index1 + index2 == eit.m_source + eit.m_target))
						cutThoseEdges_sampleGraph.insert(*e);
				}
			}

			boost::property_map < MyWeirdGraph, boost::edge_weight_t >::type w = get(boost::edge_weight, *graphs[seedIdx][v]);

			auto wp = weights.begin();
			for (boost::tie(e, e_end) = boost::edges(*graphs[seedIdx][v]); e != e_end; ++e)
				w[*e] = *wp++;

			int myIdx = clusterAndPtToIndex[seedIdx][std::make_pair(v, g[v].clusterPoints.size() / 2)];
			int tag = myIdx; //if it's a loop, it's adjusted
							 //special case: if it's a loop
			if (seedType[seedIdx] == ST_LOOP)
			{
				for (auto edge : cutThoseEdges_sampleGraph)
				{
					w[edge] = std::numeric_limits<double>::max();
				}
				tag = myIdx + 100000 * (j + 1);
			}

			parent[seedIdx][tag] = std::vector<vertex_descriptor>(N);
			D[seedIdx][tag] = std::vector<double>(N, std::numeric_limits<double>::max());
			auto indexMap = boost::get(boost::vertex_index, *graphs[seedIdx][v]);
			auto predMap = make_iterator_property_map(parent[seedIdx][tag].begin(), indexMap);
			auto distMap = make_iterator_property_map(D[seedIdx][tag].begin(), indexMap);
			boost::dijkstra_shortest_paths(*graphs[seedIdx][v], myIdx, predecessor_map(predMap).distance_map(distMap));
		}
	}
	std::cout << "done in " << double(begin + clock()) / CLOCKS_PER_SEC << " seconds. " << std::endl;


}

bool Distances::reconstructPath(int fixedVertex, int fixedVertexSample, int embeddedVertex, int embeddedSample, int seedIdx, int chainIdx, std::vector<std::pair<int, int>>& path)
{
	if (parent[seedIdx].empty())
		return false;

	int v1 = fixedVertex, v2 = embeddedVertex, clP1 = fixedVertexSample, clP2 = embeddedSample;

	int tag = clusterAndPtToIndex[seedIdx][std::make_pair(v1, clP1)];
	if (seedType[seedIdx] == ST_LOOP)
		tag += 100000 * (chainIdx + 1);

	bool inverse = D[seedIdx].find(tag) == D[seedIdx].end();

	if (inverse)
	{
		std::swap(v1, v2);
		std::swap(clP1, clP2);
	}

	int startVtx = clusterAndPtToIndex[seedIdx][std::make_pair(v1, clP1)];
	int startVtxTag = startVtx;
	if (seedType[seedIdx] == ST_LOOP)
		startVtxTag += 100000 * (chainIdx + 1);

	int endVtx = clusterAndPtToIndex[seedIdx][std::make_pair(v2, clP2)];

	int curVtx = endVtx;
	std::vector<std::pair<int, int>> result;

	while ((curVtx != startVtx) && (result.size() < 10000))
	{
		auto clusterAndPt = indexToClusterAndPt[seedIdx][curVtx];
		result.push_back(clusterAndPt);
		int nextCandidate = parent[seedIdx][startVtxTag][curVtx];
		curVtx = nextCandidate;
	}
	result.push_back({ v1,clP1 });
	path = result;

	if (!inverse)
		std::reverse(path.begin(), path.end());

	if ((result.size() >= 10000) && (result[0] == result[1]))
	{
		std::cout << "CRASHED: " << fixedVertex << ", " << embeddedVertex << std::endl;
		return false;
	}
	return true;
}

MyPolyline Distances::convertToActualPolyline(const G& g, const std::vector<std::pair<int, int>>& vec, const std::vector<MyPolyline>& polys, std::vector<double>& radii)
{
	MyPolyline result;
	if (vec.empty())
		return result;

	for (int i = 0; i + 1 < vec.size(); ++i)
	{
		MyPolyline chunk;
		auto distAndInfo = distanceBetweenSamples(vec[i].first, vec[i].second, vec[i + 1].first, vec[i + 1].second);

		if (distAndInfo.second.first.curve == distAndInfo.second.second.curve)
		{
			int curve = distAndInfo.second.first.curve;

			double idxStart = distAndInfo.second.first.segmentIdx, idxEnd = distAndInfo.second.second.segmentIdx;

			//std::cout << "v" << vec[i].first << " v" << vec[i + 1].first << ", curve = " << curve << ", segments = " << idxStart << ", " << idxEnd << " (dist = " << distAndInfo.first << ")"<< std::endl;

			if (std::ceil(std::min(idxStart, idxEnd)) < std::floor(std::max(idxStart, idxEnd)))
				chunk.insert(chunk.end(), polys[curve].begin() + std::ceil(std::min(idxStart, idxEnd)), polys[curve].begin() + std::floor(std::max(idxStart, idxEnd)));
			if (idxStart > idxEnd)
				std::reverse(chunk.begin(), chunk.end());
			chunk.insert(chunk.begin(), g[vec[i].first].clusterPoints[vec[i].second].p);
			result.insert(result.end(), chunk.begin(), chunk.end());

			for (int j = 0; j < chunk.size(); ++j)
			{
				double t = j / chunk.size();
				double r = g[vec[i].first].width*(1 - t) + g[vec[i + 1].first].width*t;
				radii.push_back(r);
			}
		}
		else
		{
			const double h = 0.1;
			Eigen::Vector2d p1 = g[vec[i].first].clusterPoints[vec[i].second].p, p2 = g[vec[i + 1].first].clusterPoints[vec[i + 1].second].p;
			double dist = (p1 - p2).norm();
			int nSamples = std::max((int)(dist / h), 1);
			for (int j = 0; j < nSamples; ++j)
			{
				double t = (double)j / nSamples;
				Eigen::Vector2d pp = p1*(1 - t) + p2*t;
				result.push_back(pp);
				radii.push_back(g[vec[i].first].width*(1 - t) + g[vec[i + 1].first].width*t);
			}

		}
	}
	result.push_back(g[vec.back().first].clusterPoints[vec.back().second].p);
	radii.push_back(g[vec.back().first].width);
	return result;
}

std::tuple<std::vector<MyPolyline>, std::vector<std::vector<double>>, std::vector<std::array<bool, 2>>, std::vector<std::array<bool, 2>>, std::vector<std::pair<PointOnCurve, PointOnCurve>> > topoGraphEmbedding(G & g, const std::vector<MyPolyline>& polys, const cv::Mat& bwImg)
{
	std::vector<std::pair<PointOnCurve, PointOnCurve>> yJunctionInfo;
	std::clock_t begin = -std::clock();
	std::cout << "[topoGraphEmbedding]: starting... ";
	std::map<edge_descriptor, size_t> myChain;
	auto chains = chainDecomposition(g, myChain);

	std::vector<std::vector<edge_descriptor>> chainsSeparated;
	auto tGraph = topoGraphHighValenceSeparated(g, chainsSeparated);
	Distances d(g, tGraph, chainsSeparated, bwImg);

	std::cout << "TOPO GRAPH: ";
	for (const auto& e : tGraph)
		std::cout << e.first << " - " << e.second << std::endl;

	std::cout << "Special vertices for deg-3 separation: ";
	std::set<size_t> separationVertices;
	for (const auto& e : tGraph)
	{
		if (boost::degree(e.first, g) == 2)
			separationVertices.insert(e.first);
		if (boost::degree(e.second, g) == 2)
			separationVertices.insert(e.second);
	}
	for (size_t s : separationVertices)
		std::cout << s << " ";
	std::cout << std::endl;

	std::map<size_t, std::pair<size_t, int>> result;
	std::vector<MyPolyline> actualResult;
	std::vector<std::vector<double>> radii;
	std::vector<std::array<bool, 2>> protectedEnds, isItASpecialDeg2Vertex;

	for (int seedIdx = 0; seedIdx < d.seeds.size(); ++seedIdx)
	{
		if (d.vertexSets[seedIdx].size() < 3) //nothing to optimize here
			continue;

		size_t theVertex = d.seeds[seedIdx];

		if (d.seedType[seedIdx] == Distances::ST_VAL2) //nothing to optimize, just find the shortest path
		{
			int chain = d.adjChains[theVertex][0];
			size_t v1 = chainsSeparated[chain].front().m_source;
			size_t v2 = chainsSeparated[chain].back().m_target;
			std::vector<std::pair<int, int>> path;
			radii.push_back({});
			bool result = d.reconstructPath(v1, g[v1].clusterPoints.size() / 2, v2, g[v2].clusterPoints.size() / 2, seedIdx, -1, path);
			actualResult.push_back(d.convertToActualPolyline(g, path, polys, radii.back()));
			protectedEnds.push_back({ false,false });
			continue;
		}

		std::vector<std::pair<size_t, int>> potentialLocations;
		std::cout << std::endl << "Starting component " << seedIdx << ", seedPt: " << theVertex << " (size: " << d.vertexSets[seedIdx].size() << "), ";

		switch (d.seedType[seedIdx])
		{
		case Distances::ST_REGULAR:
			std::cout << "a regular seed" << std::endl;
			break;
		case Distances::ST_LOOP:
			std::cout << "a seed with a loop" << std::endl;
			break;
		case Distances::ST_VAL2:
			std::cout << "a valence 2 seed" << std::endl;
			break;
		case Distances::ST_ADJ_TO_ANOTHER_SEED:
			std::cout << "adjacent to another seed" << std::endl;
			break;
		}

		std::cout << std::endl << std::endl;

		for (size_t v : d.vertexSets[seedIdx])
		{
			potentialLocations.push_back(std::make_pair(v, g[v].clusterPoints.size() / 2)); //todo: replace with a median
		}

		std::cout << potentialLocations.size() << " locs...";

		std::vector<std::pair<size_t, size_t>> adjTopoEdges;
		std::vector<int> myAdjChains;

		for (int j = 0; j < tGraph.size(); ++j)
		{
			auto& e = tGraph[j];
			if ((e.first == theVertex) || (e.second == theVertex))
			{
				if (boost::degree(e.second, g) > 2)
					std::swap(e.first, e.second);

				adjTopoEdges.push_back(e);
				myAdjChains.push_back(j);
			}
		}

		if ((d.seedType[seedIdx] == Distances::ST_REGULAR) || (d.seedType[seedIdx] == Distances::ST_LOOP))
		{
			double test = 0;

			std::vector<double> sums(potentialLocations.size(), 0);
			for (int ii = 0; ii < adjTopoEdges.size(); ++ii)
			{
				auto e = adjTopoEdges[ii];
				std::cout << e.first << "-" << e.second << std::endl;
				size_t tag = d.clusterAndPtToIndex[seedIdx][std::make_pair(e.second, g[e.second].clusterPoints.size() / 2)];

				if (d.seedType[seedIdx] == Distances::ST_LOOP)
					tag += 100000 * (myAdjChains[ii] + 1);

				for (int i = 0; i < potentialLocations.size(); ++i)
				{
					int idx = d.clusterAndPtToIndex[seedIdx][potentialLocations[i]];
					double w = d.D[seedIdx][tag][idx];

					std::vector<size_t> endPts;
					bool closeToAnEndPoint = false;
					for (const auto& ee : adjTopoEdges)
					{
						for (size_t vv : {ee.first, ee.second})
						{
							if ((boost::degree(vv, g) == 1) && (potentialLocations[i].first != vv) && ((g[vv].location - g[potentialLocations[i].first].location).norm() < 1))
								closeToAnEndPoint = true;
						}
					}

					if (closeToAnEndPoint)
					{
						w = 1000;
					}

					if (w > 1e5)
						w = 1000;

					if ((potentialLocations[i].first == e.second) && (d.seedType[seedIdx] == Distances::ST_LOOP)) //corner case: loops minimum of distance is when they are collapsed
						w = 1000;

					sums[i] += w;
				}
			}

			double myMin = 1000;
			int myWinner = -1;
			for (int i = 0; i < potentialLocations.size(); ++i)
			{
				if (myMin > sums[i])
				{
					myMin = sums[i];
					myWinner = i;
				}
			}

			result[theVertex] = potentialLocations[myWinner];
			std::cout << "Vertex " << theVertex << " decided to attach to " << potentialLocations[myWinner].first << " vertex " << potentialLocations[myWinner].second << " sample" << std::endl;
		}
		else
		{
			//if it's a non-standard seed, the embedding is fixed
			for (const auto& e : tGraph)
			{
				if ((e.first == theVertex) || (e.second == theVertex))
				{
					result[e.first] = { e.first,g[e.first].clusterPoints.size() / 2 };
					result[e.second] = { e.second, g[e.second].clusterPoints.size() / 2 };
				}
			}
		}

		std::cout << "Reconstructing the path... ";
		std::vector<std::vector<std::pair<int, int>>> paths;

		std::vector<MyPolyline> newPolylines;
		std::vector < std::vector<double>> newRadii;
		for (int i = 0; i < adjTopoEdges.size(); ++i)
		{
			auto e = adjTopoEdges[i];
			int j = myAdjChains[i]; //matters only if it's a loop
			std::vector<std::pair<int, int>> path;
			//e.first = theVertex				
			d.reconstructPath(e.second, g[e.second].clusterPoints.size() / 2, result[e.first].first, result[e.first].second, seedIdx, j, path);
			paths.push_back(path);
			newRadii.push_back({});
			newPolylines.push_back(d.convertToActualPolyline(g, path, polys, newRadii.back()));

			std::cout << newPolylines.back().size() << std::endl;

		}

		//also, merge adjChains[0] and its continuation
		double bestDot = 1;
		int connect0ChainWith = -1;
		const int step = 10;
		Eigen::Vector2d v0 = newPolylines[0][std::max(0, (int)newPolylines[0].size() - step)] - newPolylines[0].back();
		v0.normalize();
		for (int i = 1; i < newPolylines.size(); ++i)
		{
			Eigen::Vector2d vi = newPolylines[i][std::max(0, (int)newPolylines[i].size() - step)] - newPolylines[i].back();
			double dot = vi.normalized().dot(v0);
			if (dot < bestDot)
			{
				bestDot = dot;
				connect0ChainWith = i;
			}
		}

		PointOnCurve junctionPt;

		if (connect0ChainWith != -1)
		{
			std::pair<int, int> bestPair;
			double bestDist = std::numeric_limits<double>::max();
			for (int end1 : { 0, (int)newPolylines[0].size() - 1 })
				for (int end2 : {0, (int)newPolylines[connect0ChainWith].size() - 1})
				{
					double dist = (newPolylines[0][end1] - newPolylines[connect0ChainWith][end2]).squaredNorm();
					if (dist < bestDist)
					{
						bestDist = dist;
						bestPair = std::make_pair(end1, end2);
					}
				}

			if (bestPair.first == 0)
			{
				std::reverse(paths[0].begin(), paths[0].end());
				std::reverse(newPolylines[0].begin(), newPolylines[0].end());
				std::reverse(newRadii[0].begin(), newRadii[0].end());
			}

			if (bestPair.second != 0)
			{
				std::reverse(paths[connect0ChainWith].begin(), paths[connect0ChainWith].end());
				std::reverse(newPolylines[connect0ChainWith].begin(), newPolylines[connect0ChainWith].end());
				std::reverse(newRadii[connect0ChainWith].begin(), newRadii[connect0ChainWith].end());
			}

			junctionPt.curve = actualResult.size();
			junctionPt.segmentIdx = newPolylines[0].size();
			newPolylines[0].insert(newPolylines[0].end(), newPolylines[connect0ChainWith].begin(), newPolylines[connect0ChainWith].end());
			protectedEnds.push_back({ boost::degree(paths[0].front().first, g) != 1, boost::degree(paths[connect0ChainWith].back().first, g) != 1 });
			newRadii[0].insert(newRadii[0].end(), newRadii[connect0ChainWith].begin(), newRadii[connect0ChainWith].end());
		}
		else
			std::cout << "SOMETHIGN WENT WRONG" << std::endl;

		for (int i = 1; i < adjTopoEdges.size(); ++i)
		{
			if (i == connect0ChainWith)
				continue;

			auto e = adjTopoEdges[i];

			bool keepOtherEnd = boost::degree(e.second, g) != 1;
			if (paths[i].back().first == e.second)
				protectedEnds.push_back({ false,keepOtherEnd });
			else
				protectedEnds.push_back({ keepOtherEnd,false });
		}

		for (int i = 0; i < newPolylines.size(); ++i)
		{
			if (i != connect0ChainWith)
			{
				PointOnCurve pt;
				pt.curve = actualResult.size();
				pt.segmentIdx = (int)newPolylines[i].size() - 1;
				actualResult.push_back(newPolylines[i]);
				radii.push_back(newRadii[i]);

				if (i != 0)
					yJunctionInfo.push_back(std::make_pair(junctionPt, pt));
			}
		}

		std::cout << "done." << std::endl;
	}

	std::cout << "All done in " << double(begin + clock()) / CLOCKS_PER_SEC << " seconds. " << std::endl;

	isItASpecialDeg2Vertex.resize(actualResult.size());
	for (size_t s : separationVertices)
	{
		auto p = g[s].clusterPoints[g[s].clusterPoints.size() / 2];
		std::vector<PointOnCurve> pts;
		for (int i = 0; i < actualResult.size(); ++i)
		{
			if ((p.p - actualResult[i][0]).norm() < 1e-6)
				isItASpecialDeg2Vertex[i][0] = true;

			if ((p.p - actualResult[i].back()).norm() < 1e-6)
				isItASpecialDeg2Vertex[i][1] = true;

		}
	}

	return std::tie(actualResult, radii, protectedEnds, isItASpecialDeg2Vertex, yJunctionInfo);
}
