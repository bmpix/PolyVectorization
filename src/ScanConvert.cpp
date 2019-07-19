#include "stdafx.h"
#include "ScanConvert.h"

std::set<std::array<int, 2>> scanConvert(const std::vector<Eigen::Vector2d>& poly, const cv::Mat & origMask, const std::array<Eigen::MatrixXcd, 2>& roots, int startIdx, std::map<std::array<int, 2>, 
	bool>& isRootMatchingOK, std::array<bool, 2>& hitSingularity, const std::set<std::array<int, 2>>& singularities, std::array<double,2> segment)
{
	std::set<std::array<int, 2>> pixels;

	int m = origMask.rows, n = origMask.cols;
	Eigen::Vector2d p0 = poly[startIdx];
	std::array<int, 2> pi0 = { (int)std::round(p0.y()), (int)std::round(p0.x()) };
	isRootMatchingOK[pi0] = true;
	hitSingularity = { false,false };
	Eigen::Vector2d r0 = _toEig(roots[0](pi0[0], pi0[1]));

	int dirIdx = 0;
	for (int dir : {-1, 1})
	{
		std::array<int, 2> prevPixel = pi0;
		std::array<int, 2> pi = pi0;
		for (int i = startIdx; (i >= segment[0]) && (i <= segment[1]); i += dir)
		{
			Eigen::Vector2d p = poly[i];
			pi = { (int)std::round(p.y()), (int)std::round(p.x()) };

			if (singularities.find(pi) != singularities.end())
			{
				hitSingularity[dirIdx] = true;
				break;
			}

			if (i == startIdx || (prevPixel != pi))
			{
				for (int i = -1; i <= 1; ++i)
				{
					for (int j = -1; j <= 1; ++j)
					{
						std::array<int, 2> newPixel = { pi[0] + i,pi[1] + j };
						if (!inBounds(newPixel, origMask) || singularities.find(newPixel) != singularities.end())
							continue;

						isRootMatchingOK[newPixel] = !(isRootMatchingOK[prevPixel] ^ rootMatching(roots, newPixel, prevPixel));

						pixels.insert(newPixel);
					}
				}
			}
			prevPixel = pi;
		}
		++dirIdx;
	}
	return pixels;

}
