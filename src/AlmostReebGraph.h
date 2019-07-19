#pragma once
#include "typedefs.h"
#include "graph_typedefs.h"
#include <set>

typedef std::map<std::array<int, 2>, std::vector<Eigen::Vector2d>> VectorsPerPixel;

G computeAlmostReebGraph(const cv::Mat & origMask, const std::array<Eigen::MatrixXcd, 2>& roots, const std::vector<MyPolyline>& polys, 
	std::map<std::array<int, 2>, std::vector<PixelInfo>>& pixelInfo, const std::set<std::array<int, 2>>& singularities, const Eigen::MatrixXi& indices, const Eigen::VectorXcd& X, const std::vector<std::array<bool, 2>>& endedWithASingularity);
Cluster createCluster(int seedCurve, int seedPtIdx, const cv::Mat & origMask, const std::array<Eigen::MatrixXcd, 2>& roots, 
	const std::vector<MyPolyline>& polys, const std::map<std::array<int, 2>, std::vector<PixelInfo>>& pixelInfo, const std::set<int>& onlyTheseCurves, const std::set<std::array<int, 2>>& singularities, const Eigen::MatrixXi& indices, const Eigen::VectorXcd& X);
void contractSingularityBranches(G& g);
void connectStuffAroundSingularities(G& g, const cv::Mat & origMask, const std::vector<MyPolyline>& polys, const std::set<std::array<int, 2>>& singularities, const std::array<Eigen::MatrixXcd, 2>& roots, const std::vector<std::array<bool, 2>>& endedWithASingularity);
void simpleThresholds(G& g);
