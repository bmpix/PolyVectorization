#include "stdafx.h"
#include "traceAuto.h"
#include "greedyTrace.h"
#include <iostream>
#include "opencv2/core/eigen.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <fstream>

std::vector<MyPolyline> traceAll(const cv::Mat & bwImg, const cv::Mat & origMask, const cv::Mat & extMask, const std::array<Eigen::MatrixXcd, 2>& roots, const Eigen::VectorXcd & X, const Eigen::MatrixXi& indices, std::map<std::array<int, 2>, std::vector<PixelInfo>>& pixelInfo, std::vector<std::array<bool,2>>& endedWithASingularity)
{
	std::vector<MyPolyline> result;
	std::vector<std::array<bool, 2>> protectedEnds;

	int m = bwImg.rows, n = bwImg.cols;
	std::map<std::array<int, 2>, MyPolyline> resultPerPixel;
	std::map<std::array<int, 2>, int> chosenRootIdx;
	int curveIdx = 0;
	//before enabling the pragma, look again at curveIdx
//#pragma omp parallel for
	for (int i = 0; i < m; ++i)
	{
		if (i%10 == 0)
			printf("%d/%d\r", i, m);

		for (int j = 0; j < n; ++j)
		{
			if ((origMask.at<uchar>(i, j) == 0) || ((abs(roots[0](i,j))<1e-10) && (abs(roots[1](i, j))<1e-10)))
				continue;

			std::map<std::array<int, 2>, std::vector<PixelInfo>> outNewPixelInfo;
			std::pair<MyPolyline, std::array<bool, 2>> tracingResult;

			std::map<std::array<int, 2>, CenterFit> emptyFits;
			auto root = std::abs(roots[0](i, j)) > std::abs(roots[1](i, j)) ? roots[0](i, j) : roots[1](i, j);
			//next center found, start tracing
			Eigen::Vector2d seedCenter((double)j, (double)i);
			std::array<double, 2> distToSingularity = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
			tracingResult = greedyTrace(origMask, roots, seedCenter, _toEig(root), X, pixelInfo, outNewPixelInfo, indices, curveIdx, false, distToSingularity);

			if (tracingResult.first.size() > 2) //otherwise it's just junk, man
			{
				//check it's not a duplicate
				std::set<int> candidateCurves;
				for (auto &it : outNewPixelInfo)
				{
					for (auto &record : it.second)
					{
						for (auto &pi : pixelInfo[it.first])
							candidateCurves.insert(pi.curve);
					}
					//candidates only for one pixel are more than enough
					break;
				}

				for (auto &it : outNewPixelInfo)
				{
					for (auto &record : it.second)
					{
						std::map<int, double> bestDist;

						for (auto &pi : pixelInfo[it.first])
						{
							if (candidateCurves.find(pi.curve) == candidateCurves.end())
								continue;

							double dist = (pi.p - record.p).squaredNorm();
							if (bestDist.find(pi.curve) != bestDist.end())
								bestDist[pi.curve] = std::min(bestDist[pi.curve], dist);
							else
								bestDist[pi.curve] = dist;
						}

						for (auto &ii : bestDist)
						{
							if (ii.second > 0.01)
								candidateCurves.erase(ii.first);
						}

						if (candidateCurves.empty())
							break;
					}

					if (candidateCurves.empty())
						break;
				}

				if (!candidateCurves.empty())
					continue;

				resultPerPixel[{i, j}] = tracingResult.first;
				endedWithASingularity.push_back({ distToSingularity[0] < 1e10, distToSingularity[1] < 1e10 });
				curveIdx++;

				for (auto& it : outNewPixelInfo)
				{
					pixelInfo[it.first].insert(pixelInfo[it.first].end(), it.second.begin(), it.second.end());
				}

			}
		}
	}

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (resultPerPixel.find({ i,j }) == resultPerPixel.end())
				continue;

			result.insert(result.end(), resultPerPixel[{i, j}]);
		}

	std::cout << "Done. " << result.size() << " curves" << std::endl;
	return result;
}