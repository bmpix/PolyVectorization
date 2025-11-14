#pragma once
#include <QWidget>
#include "ui_mainwindow.h"
#include <QtWidgets>
#include <array>
#include <Eigen/Dense>
#include "imagequiverviewer.hpp"
#include "typedefs.h"
#include "MyScrollArea.h"
class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	MainWindow();
	~MainWindow();
	void setImage(QString filename, const cv::Mat& mask);
	void setRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);
	void setExtraRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);
	void setPolys(const std::vector<MyPolyline>& newPolys);
	void setGraph(std::string s, const G & g);
	void setSplitVertices(const std::set<size_t>& vtx);
	void setIncontractibleLoops(const std::vector<std::vector<cv::Point2f>>& loops);
	void setVectorization(std::string s, const std::vector<MyPolyline>& vectorization);
protected:
	void keyPressEvent(QKeyEvent * event);
	void redraw();
	void setScale(double newScale);
	virtual void wheelEvent(QWheelEvent * event);
private:
	const int numGraphs, numVectorizations;
	Ui::MainWindow ui;
	QImage image;
	ImageQuiverViewer *imageLabel;
	MyScrollArea *scrollArea;
	double scale;
	int origWidth, origHeight;
	std::vector<QCheckBox*> checkGraphs, checkVectorizations;
	QCheckBox* checkPolys, *checkIncontractibleLoops;
	QRadioButton* maskRadio, *imageRadio;
	QImage actualImage, maskImage;
	std::map<std::string, int> graphNameToIdx, vectorizationNameToIdx;
	QVBoxLayout *l;
};
