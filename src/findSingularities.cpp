#include "stdafx.h"
#include "findSingularities.h"
#include <iostream>
#include <set>

bool isItASingularity(const Eigen::Vector2d& p, const std::array<std::complex<double>,2> myRoots, const std::array<Eigen::MatrixXcd, 2>& roots, const cv::Mat& origMask)
{
	//return false;
	std::vector<std::array<int, 2>> shifts = { { -1,-1 },{ -1,0 },{ 0,-1 },{ 0,0 } };
	int i = std::round(p.y()), j = std::round(p.x());
	for (auto shift : shifts)
	{
		std::vector<std::array<std::complex<double>, 2>> pixelRoots;
		std::vector<std::array<int, 2>> pixelOrder = { { 0,0 },{ 0,1 },{ 1,1 },{ 1,0 } };
		for (auto v : pixelOrder)
		{
			int i1 = i + shift[0] + v[0], j1 = j + shift[1] + v[1];
			if (inBounds(i1, j1, origMask))
			{
				if ((i == i1) && (j == j1))
					pixelRoots.push_back(myRoots);
				else
					pixelRoots.push_back({ roots[0](i1,j1),roots[1](i1,j1) });
			}
		}
		if (pixelRoots.size() == 4)
		{
			bool overallMatch = true;
			for (int k = 0; k < 4; ++k)
			{
				bool match = rootMatching(pixelRoots[k], pixelRoots[(k + 1) % 4]);
				overallMatch = !(match ^ overallMatch);
			}
			if (overallMatch == false)
			{
				return true;

			}
		}
	}
	return false;
}

std::set<std::array<int,2>> findSingularities(std::array<Eigen::MatrixXcd, 2>& roots, Eigen::VectorXcd& X, const Eigen::MatrixXi& indices, const cv::Mat& origMask)
{
	std::cout << "Singularities: ";
	std::set<std::array<int, 2>> singularities;
	int m = origMask.rows, n = origMask.cols;
	for (int i=0; i<m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (isItASingularity({ j,i }, { roots[0](i,j),roots[1](i,j) }, roots, origMask))
			{
				if (singularities.find({ i,j }) == singularities.end())
					std::cout << i << ", " << j << "; ";
				singularities.insert({ i,j });
			}
		}
	std::cout << std::endl;
	return singularities;
}