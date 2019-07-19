#ifndef _GRAPH_TYPEDEFS_H_
#define _GRAPH_TYPEDEFS_H_

#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include "simple_svg_1.0.0.hpp"
#include "typedefs.h"

struct Cluster
{
	Eigen::Vector2d location;
	mutable Eigen::Vector2d root;
	std::vector<PointOnCurve> clusterPoints;
	std::vector<Eigen::Vector2d> clusterCurve;
	int clusterIdx;
	int seedCurve;
	bool sharpCorner; //filled in SplitEmUp
	bool clusterCurveHitSingularity;
	bool nextToSingularity; //if one of the cluster points is right next to a curve end which is at singularity
	double width;
	bool split;
};

struct Edge
{
	int edgeCurve;
	double weight;
};

//typedef boost::property<boost::vertex_index_t, int> VertexProperty;
using G = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Cluster, Edge>;
using listG = boost::adjacency_list<boost::vecS, boost::listS, boost::undirectedS, Cluster, Edge>;
using V = G::vertex_descriptor;
using E = std::pair<V, V>;
using vertex_descriptor = boost::graph_traits<G>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<G>::edge_descriptor;
using edge_iter = boost::graph_traits<G>::edge_iterator;
using oedge_iter = boost::graph_traits<G>::out_edge_iterator;

std::set<edge_descriptor> contract_edges(const std::vector<std::vector<boost::graph_traits<G>::edge_descriptor>>& edgeLoops, G& g);

G listToVec(const listG& g);

#endif