#pragma once
#include <QtWidgets>
#include <Eigen/Dense>
#include <array>
#include "typedefs.h"
#include "graph_typedefs.h"
#include <vector>
#include "TopoGraphEmbedding.h"
class ImageQuiverViewer : public QLabel {
	Q_OBJECT

public:
	ImageQuiverViewer(QWidget * parent = Q_NULLPTR);
	~ImageQuiverViewer();

	template <typename T>
	QPointF modelToScreen(const T & p)
	{
		T result = (p + T(0.5, 0.5))*scale;
		return{ result.x(), result.y() };
	}
	void setRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);
	void setExtraRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);
	void setPolys(const std::vector<MyPolyline>& newPolys);
	void drawRoots(const std::array<Eigen::MatrixXcd, 2>& thoseOnes, QPainter & painter);
	void drawGraphs(QPainter& painter);
	void drawPolys(QPainter & painter, const std::vector<MyPolyline>& curves, double width, QColor color);
public:
	double scale;
	std::vector<bool> graphHidden;
	std::vector<G> graphs;
	std::vector<std::string> graphTags;
	std::vector<bool> vectorizationsHidden;
	bool polysHidden;
	bool showCircles;
	bool showIncontractibleLoops;
	
	std::set<size_t> splitVtx;
	std::vector<std::vector<QPointF>> incontractibleLoops;
	std::unique_ptr<Distances> tmpDistances;
	std::vector<std::vector<MyPolyline>> vectorizations;
protected:
	bool areGraphsVisible();
	void paintEvent(QPaintEvent * event);
	void drawClusters(QPainter & painter);
	virtual void mousePressEvent(QMouseEvent * event);
private:
	std::array<Eigen::MatrixXcd, 2> roots, extraRoots;
	std::vector<MyPolyline> polys;
	std::set<std::pair<int, int>> clustersToDraw;
	std::vector<QColor> colors;
};
