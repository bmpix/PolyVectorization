#include "stdafx.h"
#include "greedyTrace.h"
#include "FindRoots.h"
#include "chooseRoot.h"
#include "bilinearInterpCoeff.h"
#include <iostream>
#include <set>
#include "findSingularities.h"

std::pair<MyPolyline, std::array<bool, 2>> greedyTrace(const cv::Mat & origMask, const std::array<Eigen::MatrixXcd, 2>& allRoots, 
	const Eigen::Vector2d & seedCenter, const Eigen::Vector2d& initialDir, const Eigen::VectorXcd & X, const std::map<std::array<int, 2>, 
	std::vector<PixelInfo>>& pixelInfo, std::map<std::array<int, 2>, std::vector<PixelInfo>>& outNewPixelInfo, const Eigen::MatrixXi& indices, 
	 int myCurveIdx, bool perpendicular, std::array<double,2>& closestDistToSingularity)
{
	//a word of advice to those modifying this code: unfortunately, this h is not a global costant, although it should be. That means that in *many* places in the code there are numbers like 10 or 0.1
	//those should actually be just h or 1/h, but I didn't fix it - sorry! so unless you really need to modify this h, don't, but if you do, just search for any constant that looks like h or 1/h and fix that too.
	const double h = 0.1;
	using namespace cv;
	Rect bigBounds(Point(), origMask.size());
	int m = origMask.rows, n = origMask.cols;
	std::array<bool, 2> protectedEnds;
	MyPolyline result;

	int firstDirectionPolySize;
	std::set<std::array<int,2>> singularities;
	for (int idx = 0; idx < 2; ++idx)
	{
		
		bool finishedWithAnotherCurve = false;
		Eigen::Vector2d p = seedCenter;
		//std::cout << p.x() << " " << p.y() << std::endl;
		std::array<int, 2> p0i = { static_cast<int>(std::round(p.y())), static_cast<int>(std::round(p.x())) };
		//std::cout << allRoots[0](p0i[0], p0i[1]) << " --- " << allRoots[1](p0i[0], p0i[1]) << std::endl;
		MyPolyline polyline;
		
		Eigen::Vector2d dir = (2*idx-1)*chooseRoot(allRoots[0](p0i[0],p0i[1]), allRoots[1](p0i[0], p0i[1]), initialDir); //sign changes depending on idx
		Eigen::Vector2d perp(dir.y(), -dir.x());
		std::vector<PixelInfo> pixelInfoQueue;

		auto emptyQueues = [&](int kk)
		{
			for (int k = 0; k < kk; ++k)
			{
				int i = std::round(pixelInfoQueue[k].p.y()), j = std::round(pixelInfoQueue[k].p.x());
				if ((i != p0i[0]) || (j != p0i[1]) || (idx == 1))
				{
					if (idx == 0)
						pixelInfoQueue[k].segmentIdx = -pixelInfoQueue[k].segmentIdx; //so that after we trace both directions, we adjust the segmentIdx correctly

					outNewPixelInfo[{i, j}].push_back(pixelInfoQueue[k]);
				}
			}
			pixelInfoQueue.erase(pixelInfoQueue.begin(), pixelInfoQueue.begin()+kk);
			//std::cout << "Emptying queues" << std::endl;
		};

		while (polyline.size() < 10000)
		{
			std::array<int, 2> pi = { static_cast<int>(std::round(p.y())),static_cast<int>(std::round(p.x()))};
			if (!inBounds(p,origMask) || (abs(X[indices(pi[0],pi[1])])<1e-10))
				break;

			if ((idx == 0) || result.empty() || (polyline.size() > 1) || ((result.back() - p).norm() > 1e-6))
			{
				polyline.push_back(p);

				PixelInfo pInfo;
				pInfo.p = p;
				pInfo.tangent = perpendicular ? perp : dir;
				pInfo.curve = myCurveIdx;
				pInfo.segmentIdx = polyline.size() - 1;
				pixelInfoQueue.push_back(pInfo);
			}


			if (pixelInfoQueue.size() > 3 / h) //if we've travelled at least two-pixel distance and accumulated some coverage.
				emptyQueues(pixelInfoQueue.size()/2);

			//check if this pixel is already covered by some curve with a similar tangent
			auto it = pixelInfo.find(pi);
			auto it2 = outNewPixelInfo.find(pi);
			
			if (it != pixelInfo.end() || it2 != outNewPixelInfo.end())
			{
				const auto& myCov = it != pixelInfo.end() ? it->second : it2->second;
				std::vector<double> distances(myCov.size()), dots(myCov.size());
				for (int i = 0; i < myCov.size(); i++)
				{
					distances[i] = (p - myCov[i].p).squaredNorm();
					dots[i] = perpendicular ? fabs(myCov[i].tangent.dot(perp)) : fabs(myCov[i].tangent.dot(dir));	
				}

				int closestVertex = -1;
				double bestDist = 1e-2;

				for (int i = 0; i < myCov.size(); ++i)
				{
					//this is a very lazy way of saying this is the same root. does not affect anything, but should be rewritten properly (via matching)
					if ((dots[i]>0.98)&&(distances[i] < bestDist))
					{
						bestDist = distances[i];
						closestVertex = i;
					}
				}

				if (closestVertex != -1)
				{
					finishedWithAnotherCurve = true;
					break;
				}
			}

			auto coeffs = bilinearInterpCoeff(X, p, origMask,indices);
			auto roots = findRoots(coeffs[0], coeffs[1]);

			Eigen::Vector2d shift = chooseRoot(roots[0], roots[1], dir);

			if ((shift.norm() < 1e-6) || isItASingularity(p, roots, allRoots, origMask))
			{
				closestDistToSingularity[idx] = std::min(closestDistToSingularity[idx],h*polyline.size());
				if (singularities.find(p0i) == singularities.end())
					singularities.insert(p0i);
				break;
			}

			shift.normalize();			
			
			dir = shift;
			Eigen::Vector2d shiftPerp(shift.y(), -shift.x());
			if (!perpendicular)
				p += h*shift;
			else
				p += h*shiftPerp;
		}

		protectedEnds[idx] = finishedWithAnotherCurve;

		if (idx == 0)
		{
			std::reverse(polyline.begin(), polyline.end());
			firstDirectionPolySize = static_cast<int>(polyline.size())-1;
		}

		result.insert(result.end(), polyline.begin(), polyline.end());
		emptyQueues(pixelInfoQueue.size());
	}

	for (auto& ppi: outNewPixelInfo)
		for (auto& pinfo : ppi.second)
			pinfo.segmentIdx += firstDirectionPolySize;

	//and complete the PixelInfo by pushing the first few elements from the first curve
	//they are somewhere in the middle though
	if (!result.empty())
	{
		for (int i = 0; i <= firstDirectionPolySize; ++i)
		{
			auto p = result[i];
			std::array<int, 2> p0i = { static_cast<int>(std::round(p.y())), static_cast<int>(std::round(p.x())) };
			PixelInfo pInfo;
			pInfo.p = p;
			if (i + 1 < result.size())
				pInfo.tangent = (result[i + 1] - result[i]).normalized();
			else if (i > 0)
				pInfo.tangent = (result[i] - result[i - 1]).normalized();

			pInfo.curve = myCurveIdx;
			pInfo.segmentIdx = i;

			outNewPixelInfo[p0i].push_back(pInfo);
		}
	}

	return std::make_pair(result, protectedEnds);
}
