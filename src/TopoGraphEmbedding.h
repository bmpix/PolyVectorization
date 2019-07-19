#pragma once
#include "graph_typedefs.h"

struct Distances
{
	enum SeedType
	{
		ST_REGULAR, ST_VAL2, ST_LOOP, ST_ADJ_TO_ANOTHER_SEED
	};

	Distances(const G& g, const std::vector<std::pair<size_t, size_t>>& topoGraph, const std::vector<std::vector<edge_descriptor>>& chainsSeparated, const cv::Mat& bwImg);
	bool reconstructPath(int fixedVertex, int fixedVertexSample, int embeddedVertex, int embeddedSample, int seedIdx, int chainIdx, std::vector<std::pair<int, int>>& path);
	std::vector<std::map<std::pair<int, int>, int>> clusterAndPtToIndex;
	std::vector<std::map<int, std::pair<int, int>>> indexToClusterAndPt;
	std::vector<std::map<int,std::vector<double>>> D;
	std::vector<std::map<int, std::vector<vertex_descriptor>>> parent;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, boost::property<boost::edge_weight_t, double, 
		boost::property<boost::edge_weight2_t, double>>> MyWeirdGraph;
	const G& g;
	std::vector<std::map<size_t,std::shared_ptr<MyWeirdGraph>>> graphs;
	std::vector<std::set<size_t>> vertexSets; //
	const std::vector<std::pair<size_t, size_t>>& topoGraph;
	const std::vector<std::vector<edge_descriptor>>& chainsSeparated;
	std::vector<size_t> seeds;
	std::vector<SeedType> seedType;
	std::map<size_t, std::vector<int>> adjChains; //specifies indices of chains adjacent to that vertex
	std::map<size_t, size_t> orient(const G& g, const std::vector<int>& chains, const std::set<size_t>& vertSet, int root, const std::set<edge_descriptor>& cutThoseEdges);
	std::pair<double, std::pair<PointOnCurve, PointOnCurve>> distanceBetweenSamples(int v1, int s1, int v2, int s2);
	MyPolyline convertToActualPolyline(const G& g, const std::vector<std::pair<int, int>>& vec, const std::vector<MyPolyline>& polys, std::vector<double>& radii);
	void fillInSeedTypes();

	std::map<std::pair<size_t, int>, std::vector<size_t>> vertexAndCurveToSamples;
	std::map<std::pair<size_t, size_t>, std::vector<int>> pairOfVerticesToSharedCurves;
};

std::tuple<std::vector<MyPolyline>, std::vector<std::vector<double>>, std::vector<std::array<bool, 2>>, std::vector<std::array<bool, 2>>, std::vector<std::pair<PointOnCurve, PointOnCurve>> > topoGraphEmbedding(G& g, const std::vector<MyPolyline>& polys, const cv::Mat& bwImg);
