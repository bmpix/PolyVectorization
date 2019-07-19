#include "stdafx.h"
#include "mainwindow.h"
#include <QtWidgets>
#include "ContractDeg2.h"
MainWindow::MainWindow() :scale(1.0),numGraphs(8),numVectorizations(2) {
	//ui.setupUi(this);
	

	imageLabel = new ImageQuiverViewer;
	
	scrollArea = new MyScrollArea;
	imageLabel->setBackgroundRole(QPalette::Base);
	imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageLabel->setScaledContents(true);

	scrollArea->setBackgroundRole(QPalette::Dark);
	scrollArea->setWidget(imageLabel);
	scrollArea->setVisible(false);
	setCentralWidget(scrollArea);
	
	l = new QVBoxLayout;
	for (int i = 0; i < numGraphs; ++i)
	{
		checkGraphs.push_back(new QCheckBox);
		checkGraphs.back()->setChecked(false);
		connect(checkGraphs[i], &QCheckBox::stateChanged, [&](int ignore) {
			for (int j = 0; j < numGraphs; ++j)
			{
				imageLabel->graphHidden[j] = !checkGraphs[j]->isChecked();
			}
			imageLabel->update();
		});
	}

	checkPolys = new QCheckBox;
	checkPolys->setChecked(true);
	checkPolys->setText("Curves");
	connect(checkPolys, &QCheckBox::stateChanged, [&]() {
		imageLabel->polysHidden = !checkPolys->isChecked();
		imageLabel->update();
	});

	checkIncontractibleLoops = new QCheckBox;
	checkIncontractibleLoops->setChecked(false);
	checkIncontractibleLoops->setText("Incontractible");
	connect(checkIncontractibleLoops, &QCheckBox::stateChanged, [&]() {
		imageLabel->showIncontractibleLoops = checkIncontractibleLoops->isChecked();
		imageLabel->update();
	});

	maskRadio = new QRadioButton;
	maskRadio->setText("Mask");
	maskRadio->setChecked(false);
	connect(maskRadio, &QRadioButton::toggled, [&]() {
		redraw();
		imageLabel->update();
	});

	imageRadio = new QRadioButton;
	imageRadio->setText("Image");
	imageRadio->setChecked(true);

	l->addWidget(maskRadio);
	l->addWidget(imageRadio);

	l->addWidget(checkPolys);
	l->addWidget(checkIncontractibleLoops);

	for (int i = 0; i < numGraphs; ++i)
		l->addWidget(checkGraphs[i]);


	for (int i = 0; i < numVectorizations; ++i)
	{
		checkVectorizations.push_back(new QCheckBox);
		checkVectorizations.back()->setChecked(false);
		connect(checkVectorizations[i], &QCheckBox::stateChanged, [&](int ignore) {
			for (int j = 0; j < numVectorizations; ++j)
			{
				imageLabel->vectorizationsHidden[j] = !checkVectorizations[j]->isChecked();
			}
			imageLabel->update();
		});
		l->addWidget(checkVectorizations[i]);
	}


	QDockWidget* dockWidget = new QDockWidget(tr(""), this);
	dockWidget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	QWidget* garbage = new QWidget();
	garbage->setLayout(l);
	dockWidget->setWidget(garbage);
	addDockWidget(Qt::RightDockWidgetArea, dockWidget);

	imageLabel->graphHidden.resize(numGraphs,true);
	imageLabel->graphs.resize(numGraphs);
	imageLabel->graphTags.resize(numGraphs);
	imageLabel->polysHidden = false;
	imageLabel->showCircles = true;
	imageLabel->showIncontractibleLoops = false;
	imageLabel->vectorizationsHidden.resize(numVectorizations, true);
	imageLabel->vectorizations.resize(numVectorizations);
}

MainWindow::~MainWindow() {
	
}

QImage Mat2QImage(const cv::Mat_<uchar> &src)
{
	QImage dest(src.cols, src.rows, QImage::Format_ARGB32);
	for (int y = 0; y < src.rows; ++y) {
		const uchar *srcrow = src[y];
		QRgb *destrow = (QRgb*)dest.scanLine(y);
		for (int x = 0; x < src.cols; ++x) {
			unsigned int color = srcrow[x];
			destrow[x] = qRgba(color, color, color, 255);
		}
	}
	return dest;
}

void MainWindow::setImage(QString filename, const cv::Mat& mask)
{
	actualImage = QImage (filename);
	imageLabel->setPixmap(QPixmap::fromImage(actualImage));
	origWidth = imageLabel->pixmap()->width();
	origHeight = imageLabel->pixmap()->height();
	scale = 1.0;
	scrollArea->setVisible(true);
	imageLabel->adjustSize();

	maskImage = Mat2QImage(mask);
}

void MainWindow::setRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots)
{
	imageLabel->setRoots(newRoots);
}

void MainWindow::setExtraRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots)
{
	imageLabel->setExtraRoots(newRoots);
}

void MainWindow::setPolys(const std::vector<MyPolyline>& newPolys)
{
	imageLabel->setPolys(newPolys);
}



void MainWindow::setGraph(std::string s, const G& g)
{
	int idx = graphNameToIdx.size();
	graphNameToIdx.insert(std::make_pair(s, idx));
	checkGraphs[idx]->setText(QString::fromStdString(s));
	imageLabel->graphs[idx] = g;
	imageLabel->graphTags[idx] = s;
}

void MainWindow::setSplitVertices(const std::set<size_t>& vtx)
{
	imageLabel->splitVtx = vtx;
}

void MainWindow::setIncontractibleLoops(const std::vector<std::vector<cv::Point2f>>& loops)
{
	imageLabel->incontractibleLoops.resize(loops.size());
	for (int i = 0; i < loops.size(); ++i)
	{
		imageLabel->incontractibleLoops[i].resize(loops[i].size() + 1);
		for (int j = 0; j < loops[i].size(); ++j)
			imageLabel->incontractibleLoops[i][j] = { loops[i][j].x,loops[i][j].y };
		imageLabel->incontractibleLoops[i][loops[i].size()] = { loops[i].front().x,loops[i].front().y };
	}
}

void MainWindow::setVectorization(std::string s, const std::vector<MyPolyline>& vectorization)
{
	int idx = vectorizationNameToIdx.size();
	vectorizationNameToIdx.insert(std::make_pair(s, idx));
	checkVectorizations[idx]->setText(QString::fromStdString(s));
	imageLabel->vectorizations[idx] = vectorization;
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_C)
	{
		imageLabel->showCircles = !imageLabel->showCircles;
		imageLabel->update();
	}
	if (event->key() == Qt::Key_A)
	{
		//autoscale
		double w = centralWidget()->width();
		double h = centralWidget()->height();
		setScale(std::min(w / origWidth,  h / origHeight));
	}
	if (event->key() == Qt::Key_Plus)
	{
		setScale(scale + (10 / 300.0)*scale);
	}
	else if (event->key() == Qt::Key_Minus)
	{
		setScale(scale - (10 / 300.0)*scale);
	}
}

void MainWindow::redraw()
{
	if (scale < 10000.0 / std::max(origHeight, origHeight))
	{
		imageLabel->setFixedSize(QWIDGETSIZE_MAX, QWIDGETSIZE_MAX);
		if (maskRadio->isChecked())
			imageLabel->setPixmap(QPixmap::fromImage(maskImage).scaled(origWidth*scale, origHeight*scale, Qt::KeepAspectRatio, Qt::FastTransformation));
		else
			imageLabel->setPixmap(QPixmap::fromImage(actualImage).scaled(origWidth*scale, origHeight*scale, Qt::KeepAspectRatio, Qt::FastTransformation));
	}
	else
	{
		imageLabel->setPixmap(QPixmap());
		imageLabel->setFixedSize(origWidth*scale, origHeight*scale);
	}
	imageLabel->resize(imageLabel->pixmap()->size());
}

void MainWindow::setScale(double newScale)
{
	scale = newScale;
	scale = std::max(0.5, scale);
	imageLabel->scale = scale;
	redraw();
}

void MainWindow::wheelEvent(QWheelEvent * event)
{
	auto km = QApplication::keyboardModifiers();
	if (km & Qt::KeyboardModifier::AltModifier && imageLabel->underMouse())
	{
		double oldScale = scale;

		auto origPos = event->pos();
		auto eventPos = imageLabel->mapFromParent(event->pos());
		Eigen::Vector2d mousePos(eventPos.x(), eventPos.y());
		mousePos = mousePos / scale - Eigen::Vector2d(0.5, 0.5);

		setScale(scale + (event->delta() / 300.0)*scale);

		scrollArea->horizontalScrollBar()->setValue(scrollArea->horizontalScrollBar()->maximum()*mousePos.x() / actualImage.width());
		scrollArea->verticalScrollBar()->setValue(scrollArea->verticalScrollBar()->maximum()*mousePos.y() / actualImage.height());

		event->accept();
	}
}