#include "stdafx.h"
#include "bilinearInterpCoeff.h"

std::array<std::complex<double>, 2> bilinearInterpCoeff(const Eigen::VectorXcd & X, const Eigen::Vector2d & p, const cv::Mat & extMask, const Eigen::MatrixXi & indices)
{
	double i = p.y();
	double j = p.x();
	int i1 = std::floor(i), i2 = std::ceil(i);
	int j1 = std::floor(j), j2 = std::ceil(j);
	if (i1 == i2)
		i2 = i1 + 1;
	if (j1 == j2)
		j2 = j1 + 1;

	int nonzeros = X.size() / 2;

	int m = extMask.rows, n = extMask.cols;
	if ((i2 == m) || (j2 == n))
		return{ X(indices((int)std::round(i),(int)std::round(j))),X(indices((int)std::round(i),(int)std::round(j)) + nonzeros) };


	std::array<std::complex<double>, 2> result;

	for (int rootIdx = 0; rootIdx < 2; rootIdx++)
	{
		std::complex<double> X_11 = extMask.at<uchar>(i1, j1) != 0 ? X(indices(i1, j1) +rootIdx*nonzeros) : 0;
		std::complex<double> X_12 = extMask.at<uchar>(i1, j2) != 0 ? X(indices(i1, j2) +rootIdx*nonzeros) : 0;
		std::complex<double> X_21 = extMask.at<uchar>(i2, j1) != 0 ? X(indices(i2, j1) +rootIdx*nonzeros) : 0;
		std::complex<double> X_22 = extMask.at<uchar>(i2, j2) != 0 ? X(indices(i2, j2) +rootIdx*nonzeros) : 0;
		Eigen::Matrix2cd M;
		M << X_11, X_12, X_21, X_22;
		result[rootIdx] = (1.0 / ((i2 - i1)*(j2 - j1)))*Eigen::Vector2d(i2 - i, i - i1).dot(M*Eigen::Vector2d(j2 - j, j - j1));
	}

	return result;
}

