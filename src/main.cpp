#include "stdafx.h"
#define WITH_GUI 1
#ifdef _WIN32
#include <QtCore/QCoreApplication>
#include "mainwindow.h"
#include <QElapsedTimer>
#else
#include <sys/time.h>
#include <unistd.h>
#endif

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include "Eigen/Dense"
#include "opencv2/core/eigen.hpp"
#include "Eigen/Core"
#include <fstream>
#include "Optimizer.h"
#include <array>
#include "FindRoots.h"
#include "traceAuto.h"
#include "chopFakeEnds.h"
#include "simple_svg_1.0.0.hpp"
#include "Simplify.h"
#include "AlmostReebGraph.h"
#include "ContractLoops.h"
#include "RemoveShortBranches.h"
#include "SplitEmUp.h"
#include "ChainDecomposition.h"
#include "polynomial_energy.h"
#include "findSingularities.h"
#include "Params.h"
#include "ContractDeg2.h"
#include "Smooth.h"
#include "ContractLoops2.h"
#include "TopoGraphEmbedding.h"
cv::Mat bwImg, origMask;
int m, n;
Eigen::MatrixXcd g,tau, tauTimesGmag;
Eigen::MatrixXd gMag,weight;
Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

void calculateGradient()
{
	using namespace cv;
	Mat grad_x, grad_y;
	const double scale = 1.0;
	const double delta = 0;

	Sobel(bwImg, grad_x, CV_32F, 1, 0, 3, scale, delta, BORDER_DEFAULT);
	Sobel(bwImg, grad_y, CV_32F, 0, 1, 3, scale, delta, BORDER_DEFAULT);

	Mat planes[] = { grad_x, grad_y };
	Mat cvG;
	merge(planes, 2, cvG);
	cv2eigen(cvG, g);

	tauTimesGmag = g*std::complex<double>(0.0, 1.0);
	gMag = tauTimesGmag.cwiseAbs();

	double maxGradMag = gMag.maxCoeff();
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (fabs(gMag(i, j) / maxGradMag) < 0.1)
			{
				gMag(i, j) = 0;
				tauTimesGmag(i, j) = 0;
			}
		}

	Eigen::MatrixXd gMagNoZeros = gMag;

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (fabs(gMag(i, j)) < 1e-10)
				gMagNoZeros(i, j) = 1;
		}

	tau = tauTimesGmag.array() / gMagNoZeros.array();
}

void calculateWeight()
{
	using namespace cv;

	Eigen::MatrixXcd eigTauTimesGmag2 = tauTimesGmag.array().pow(2);

	//calculate laplacian
	//todo: get rid of so many copies back and forth
	Mat eigTauTimesGmag2Re, eigTauTimesGmag2Im;

	Eigen::MatrixXd eigTauTimesGmag2Real = eigTauTimesGmag2.real(), eigTauTimesGmag2Imag = eigTauTimesGmag2.imag();
	eigen2cv(eigTauTimesGmag2Real, eigTauTimesGmag2Re);
	eigen2cv(eigTauTimesGmag2Imag, eigTauTimesGmag2Im);
	
	Mat cv_gMagPow2;
	Eigen::MatrixXd gMagPow2 = gMag.array().pow(2);
	eigen2cv(gMagPow2, cv_gMagPow2);

	Mat Lx, Ly;
	Mat kernel;
	kernel = Mat::ones(3, 3, CV_64F);
	kernel.at<double>(1, 1) = 0;
	filter2D(eigTauTimesGmag2Re, Lx, -1, kernel);
	filter2D(eigTauTimesGmag2Im, Ly, -1, kernel);

	Eigen::MatrixXd Lx_eig, Ly_eig;
	cv2eigen(Lx, Lx_eig);
	cv2eigen(Ly, Ly_eig);

	Eigen::MatrixXcd mse = Lx_eig + std::complex<double>(0, 1)*Ly_eig;
	Eigen::MatrixXd mseNorm = mse.cwiseAbs();
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			if (mseNorm(i,j)<1e-10)
				mseNorm(i, j) = 1;

	mse = mse.array() / mseNorm.array();

	Eigen::MatrixXcd tau2 = tau.array().pow(2);
	mse = mse - tau2;

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			if (fabs(gMag(i, j)) < 1e-10)
				mse(i, j) = 0;

	weight = mse.cwiseAbs();

	weight = Eigen::MatrixXd::Ones(m, n) - weight / weight.maxCoeff();

	for (int i = 0; i<m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (fabs(gMag(i, j)) < 1e-10)
			{
				weight(i, j) = 0;
			}
		}
}

Eigen::MatrixXi calculateIndices(const cv::Mat& mask, int& nnz)
{
	Eigen::MatrixXi indices = Eigen::MatrixXi(m, n);
	indices.setConstant(-1);
	int idx = 0;
	int tt = mask.type();
	for (int j = 0; j < n; ++j)
		for (int i = 0; i<m; ++i)
		{
			if (mask.at<uchar>(i, j) != 0)
			{
				indices(i, j) = idx;
				idx++;
			}
		}
	nnz = idx;
	return indices;
}

void computeAllGeodesicDistances(const cv::Mat& mask, const Eigen::MatrixXi& indices, int nnz)
{
	std::cout << "Computing all geodesic distances... ";
	using namespace boost;
	typedef adjacency_list<vecS, vecS, undirectedS, no_property,
		property< edge_weight_t, int, property< edge_weight2_t, int > > > Graph;
	typedef std::pair < int, int >Edge;

	int V = nnz;
	std::vector<Edge> myEdges;
	int m = mask.rows, n = mask.cols;
	for (int i=0; i<m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (mask.at<uchar>(i, j) != 0)
			{
				for (int sign : {-1, 1}) //leftRight
				{
					for (int dir = 0; dir < 2; dir++) //horizontal or vertical
					{
						std::pair<int, int> neigh;
						if (useNeighbor(i, j, m, n, (dir == 1), sign == -1, mask, neigh))
						{
							int idx1 = indices(i, j);
							int idx2 = indices(neigh.first, neigh.second);
							if (idx1 > idx2)
								myEdges.push_back({idx1,idx2});
						}
					}
				}
			}
		}

	Graph g(&myEdges[0], &myEdges[0]+ myEdges.size(), V);

	property_map < Graph, edge_weight_t >::type w = get(edge_weight, g);
	graph_traits < Graph >::edge_iterator e, e_end;
	
	for (boost::tie(e, e_end) = edges(g); e != e_end; ++e)
		w[*e] = 1;

	std::vector < int >d(V, (std::numeric_limits < int >::max)());
	std::vector<std::vector<int>> D(V, std::vector<int>(V));
	johnson_all_pairs_shortest_paths(g, D, distance_map(&d[0]));
	std::cout << "done." << std::endl;
}

void repairMask(cv::Mat& origMask)
{
	std::vector<std::pair<int, int>> newGuys;
	for (int i = 0; i<m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (origMask.at<uchar>(i, j) == 0)
			{
				int nn = 0;
				for (int i1 = -1; i1 < 2; ++i1)
					for (int j1 = -1; j1 < 2; ++j1)
					{
						if ((i1 == j1) || (i1 + i < 0) || (i1 + i >= m) || (j1 + j < 0) || (j1 + j >= n))
							continue;

						if (origMask.at<uchar>(i1 + i, j1 + j) != 0)
							nn++;
					}
				if (nn >= 5)
				{
					newGuys.push_back({ i, j });
				}
			}
		}

	for (auto p : newGuys)
		origMask.at<uchar>(p.first, p.second) = 255;
}

std::string type2str(int type) {
	std::string r;

	uchar depth = type & CV_MAT_DEPTH_MASK;
	uchar chans = 1 + (type >> CV_CN_SHIFT);

	switch (depth) {
	case CV_8U:  r = "8U"; break;
	case CV_8S:  r = "8S"; break;
	case CV_16U: r = "16U"; break;
	case CV_16S: r = "16S"; break;
	case CV_32S: r = "32S"; break;
	case CV_32F: r = "32F"; break;
	case CV_64F: r = "64F"; break;
	default:     r = "User"; break;
	}

	r += "C";
	r += (chans + '0');

	return r;
}

int main(int argc, char *argv[])
{
#ifdef _WIN32
	QApplication a(argc, argv);
#endif
	using namespace cv;
	Mat image;

	if (argc != 2)
	{
		std::cout << "Usage: " << std::endl;
		std::cout << "polyvector_thing.exe FILE" << std::endl;
		return 0;
	}

	std::string filename = argv[1];
	
	std::cout << "Loading " << filename << "... ";
	image = imread(filename, IMREAD_COLOR);

	if (!image.data) // Check for invalid input
	{
		std::cout << "Could not open or find the image" << std::endl;
		return -1;
	}
	std::cout << "ok" << std::endl;

	cvtColor(image, bwImg, CV_BGR2GRAY);
	bwImg = Scalar(255) - bwImg;
	m = bwImg.rows; n = bwImg.cols; //m: height, n: width

#ifdef _WIN32
	QElapsedTimer timer;
	timer.start();
#else
	struct timeval start, end;
	long mtime, seconds, useconds;
	gettimeofday(&start, NULL);
#endif

	//fill in the mask
	threshold(bwImg, origMask, BACKGROUND_FOREGROUND_THRESHOLD, 255, THRESH_BINARY);

	for (int i=0; i<3; ++i)
		repairMask(origMask);

#ifdef WITH_GUI
	MainWindow mw;
	mw.setImage(QString::fromStdString(filename), origMask);
	mw.show();
#endif

	calculateGradient();
	calculateWeight();
	int nnz = 0;
	auto indices = calculateIndices(origMask, nnz);

	std::vector<double> ws;
	double maxRes = 0;

	std::vector<std::vector<Eigen::Vector2d>> centersForI(m), axiForI(m);
	std::vector<std::vector<double>> resForI(m);
	std::vector<std::vector<int>> notNullsForI(m);
	std::map<std::array<int, 2>, CenterFit> fits;

	std::cout << "Optimizing...";
	Eigen::VectorXcd X = optimize(bwImg, weight, tau, FRAME_FIELD_SMOOTHNESS_WEIGHT, origMask, indices);
	
	std::cout << "done. " << std::endl;

	std::cout << "Finding roots.. ";
	auto roots = findRoots(X, origMask);

	auto singularities = findSingularities(roots, X, indices, origMask);
	bool improved;
	int totalNSingularities = 0;
	do
	{
		int origSingularityCount = singularities.size();

		bool somethingNew = false;
		for (auto s : singularities)
		{
			if (weight(s[0], s[1]) > 1e-5)
			{
				somethingNew = true;
				weight(s[0], s[1]) = 0;
				totalNSingularities++;
			}
		}
		if (!somethingNew)
			break;

		X = optimize(bwImg, weight, tau, FRAME_FIELD_SMOOTHNESS_WEIGHT, origMask, indices);
		roots = findRoots(X, origMask);
		singularities = findSingularities(roots, X, indices, origMask);

		std::cout << "done (" << origSingularityCount - singularities.size() << " singularities removed" << std::endl;
		improved = origSingularityCount - singularities.size() > 0;
	} while (improved);

	
	std::cout << "Tracing... ";
	std::map<std::array<int, 2>, std::vector<PixelInfo>> pixelInfo;
	std::vector<std::array<bool, 2>> endedWithASingularity;
	auto polys = traceAll(bwImg, origMask, origMask, roots, X, indices,pixelInfo, endedWithASingularity);

	std::vector<svg::Color> colors = { { 0,0,0 }, {232,26,75}, {255,198,11}, {57,181,74}, {0,131,202}, {134,50,140}, {241,90,34}, {255,242,0}, {0,166,81}, {46,49,146}, {96,45,145}, {237,28,36}, {0,168,222}, {51,51,145}, {156,41,140}, {233,24,126}, {235,44,59}, {239,128,49}, {253,233,43}, {145,190,74}, {0,160,85}
	};

	auto reebGraph = computeAlmostReebGraph(origMask, roots, polys, pixelInfo, singularities, indices, X, endedWithASingularity);
	
#ifdef WITH_GUI
	mw.setRoots(roots);
	mw.setPolys(polys);
	mw.setGraph("Orig graph", reebGraph);
#endif

	contractSingularityBranches(reebGraph);
	simpleThresholds(reebGraph);

#ifdef WITH_GUI
	mw.setGraph("Sing contracted", reebGraph);
#endif

	connectStuffAroundSingularities(reebGraph, origMask, polys, singularities,roots,endedWithASingularity);

#ifdef WITH_GUI
	mw.setGraph("Sing verts", reebGraph);
#endif

	edge_iter eit, eend;
	for (std::tie(eit, eend) = boost::edges(reebGraph); eit != eend; ++eit)
	{
		reebGraph[*eit].weight = 1;
	}
	contractLoops(reebGraph, origMask,polys);

#ifdef WITH_GUI
	mw.setGraph("Contracted", reebGraph);
#endif

	std::map<edge_descriptor, size_t> ignore;
	removeBranchesFilter1(reebGraph,false,ignore);

#ifdef WITH_GUI
	mw.setGraph("Cut", reebGraph);
#endif
	
	splitEmUpCorrectly(reebGraph);

	bool chopped = true;
	while (chopped)
	{
		chopped = false;
		for (size_t v = 0; v < boost::num_vertices(reebGraph); ++v)
		{
			if ((boost::degree(v, reebGraph) == 2) && (reebGraph[v].nextToSingularity || reebGraph[v].clusterCurveHitSingularity))
			{
				std::vector<size_t> verts;
				oedge_iter eit, eend;
				for (std::tie(eit, eend) = boost::out_edges(v, reebGraph); eit != eend; ++eit)
				{
					verts.push_back(eit->m_target);
				}
				boost::clear_vertex(v, reebGraph);
				auto e = boost::add_edge(verts[0], verts[1], reebGraph);
				reebGraph[e.first].edgeCurve = -1;
			}
			if (chopped)
				break;
		}
	}

	chopped = true;
	while (chopped)
	{
		chopped = false;
		edge_iter eit, eend;
		for (std::tie(eit, eend) = boost::edges(reebGraph); eit != eend; ++eit)
		{
			if (boost::degree(eit->m_source, reebGraph) > 2 && boost::degree(eit->m_target, reebGraph) > 2)
			{
				boost::remove_edge(*eit, reebGraph);
				chopped = true;
				break;
			}
		}
	}

#ifdef WITH_GUI
	mw.setGraph("Extra chop", reebGraph);
#endif
	
	std::vector<MyPolyline> vectorization;
	std::vector<std::vector<double>> radii;
	std::vector<std::array<bool, 2>> protectedEnds;
	std::vector<std::pair<PointOnCurve, PointOnCurve>> yJunctions; //stores info about Y-junctions (where which curve connects)
	std::vector<std::array<bool, 2>> isItASpecialDeg2Vertex;
	std::tie(vectorization, radii, protectedEnds, isItASpecialDeg2Vertex, yJunctions) = topoGraphEmbedding(reebGraph, polys, bwImg);

#ifdef WITH_GUI
	mw.setVectorization("Optimized",vectorization);
#endif

	G wG;
	std::tie(vectorization, wG) = chopFakeEnds(vectorization, radii, protectedEnds, isItASpecialDeg2Vertex, yJunctions);

#ifdef WITH_GUI
	mw.setGraph("Lala", wG);
#endif
	
	for (std::tie(eit, eend) = boost::edges(wG); eit != eend; ++eit)
	{
		wG[*eit].weight = 1.0;//(wG[eit->m_source].location- wG[eit->m_target].location).norm();
		//std::cout << "edge between " << eit->m_source << " and " << eit->m_target << ": w = " << wG[*eit].weight << std::endl;
	}

	std::cout << "Finding cycles: ";
	std::vector<edge_descriptor> removedEdges;
	if (boost::num_edges(wG) < 350)
	{
		std::cout << "Using Tarjan's algorithm " << std::endl;
		removedEdges = contractLoops2(wG, origMask, vectorization);
	}
	else
	{
		std::cout << "Using min spanning trees algorithm " << std::endl;
		removedEdges = contractLoops(wG, origMask, vectorization);
	}

#ifdef WITH_GUI
	mw.setGraph("Lala Contracted", wG);
#endif
	
	std::vector<std::vector<std::pair<double, double>>> cutThosePieces(vectorization.size());
	for (auto e: removedEdges)
	{
		int curve = wG[e.m_source].clusterPoints[0].curve;
		if (curve == wG[e.m_target].clusterPoints[0].curve)
		{
			double s1 = wG[e.m_source].clusterPoints[0].segmentIdx, s2 = wG[e.m_target].clusterPoints[0].segmentIdx;
			cutThosePieces[curve].push_back(std::minmax(s1, s2));
		}
	}

	decltype(vectorization) newVectorization;

	for (int i = 0; i < vectorization.size(); ++i)
	{
		if (vectorization[i].empty())
			continue;

		std::vector<std::pair<double,double>> segments;
		std::sort(cutThosePieces[i].begin(), cutThosePieces[i].end(), [](const std::pair<double, double>& p1, const std::pair<double, double>& p2) {return p1.first < p2.first; });
		segments.push_back({ 0.0,0.0 });
		for (int j = 0; j < cutThosePieces[i].size(); ++j)
		{
			segments.back().second = cutThosePieces[i][j].first;
			segments.push_back({ cutThosePieces[i][j].second,0.0 });
		}
		segments.back().second = vectorization[i].size() - 1;

		for (int j = 0; j < segments.size(); ++j)
		{
			MyPolyline newPoly;
			if (fabs(segments[j].second - segments[j].first) > 1e-5)
			{
				for (int k = segments[j].first; k <= segments[j].second; ++k)
				{
					if (newPoly.empty()||(newPoly.back()-vectorization[i][k]).norm()>1e-6)
						newPoly.push_back(vectorization[i][k]);
				}
			}
			if (!newPoly.empty())
				newVectorization.push_back(newPoly);
		}
	}

	for (int i = 0; i<newVectorization.size(); ++i)
		newVectorization[i] = simplify(newVectorization[i], 1e-2);

	smooth(newVectorization);

#ifdef WITH_GUI
	mw.setVectorization("Final", newVectorization);
#endif
#ifdef _WIN32
	std::cout << "Total time: " << timer.elapsed()/1000 << " s" << std::endl;
#endif
	svg::Image bgImg(filename, n, m, -0.5, 0, 0.6);
	double minX = std::numeric_limits<double>::max(), minY = minX, maxX = std::numeric_limits<double>::min(), maxY = maxX;
	for (int i = 0; i<polys.size(); ++i)
		for (int j = 0; j < polys[i].size(); ++j)
		{
			minX = std::min(polys[i][j].x(), minX);
			minY = std::min(polys[i][j].y(), minY);
			maxX = std::max(polys[i][j].x(), maxX);
			maxY = std::max(polys[i][j].y(), maxY);
		}

	svg::Document doc(filename + ".svg", svg::Layout(svg::Dimensions(n, m), svg::Layout::TopLeft));
	doc << bgImg;
	for (int i = 0; i < newVectorization.size(); ++i)
	{
		if (!newVectorization[i].empty())
		{
			svg::Polyline poly(svg::Fill(), svg::Stroke(1, colors[1]));
			for (auto& p : newVectorization[i])
				poly << svg::Point(p.x() + 0.5, p.y() + 0.5);
			doc << poly;
		}
	}
	doc.save();

#ifdef WITH_GUI
	return a.exec();
#else
	return 0;
#endif

}
