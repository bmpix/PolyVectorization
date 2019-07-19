#include "stdafx.h"
#include "FindRoots.h"

std::array<std::complex<double>, 2> findRoots(std::complex<double> c0, std::complex<double> c2)
{
	//roots of z^4 + c2*z^2 + c0 = 0
	return{ sqrt((-c2 + sqrt(c2*c2 - 4.0 * c0)) / 2.0), sqrt((-c2 - sqrt(c2*c2 - 4.0 * c0)) / 2.0) };
}

std::array<Eigen::MatrixXcd, 2> findRoots(const Eigen::VectorXcd & X, const cv::Mat & mask)
{
	typedef std::complex<double> cmplx;

	int m = mask.rows, n = mask.cols;
	Eigen::MatrixXcd root1(m,n), root2(m,n);
	root1.setZero(); root2.setZero();


	int nonzeros = X.size() / 2;

	int idx = 0;
	for (int j=0; j<n; ++j)
		for (int i = 0; i < m; ++i)
		{
			if (mask.at<uchar>(i, j) == 0)
				continue;

			cmplx c0 = X(idx), c2 = X(idx + nonzeros);
			auto roots = findRoots(c0, c2);
			root1(i, j) = roots[0];
			root2(i, j) = roots[1];

			idx++;
		}
	return{ root1,root2 };
}
