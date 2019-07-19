#include "stdafx.h"
#include "Smooth.h"

void smooth(std::vector<MyPolyline>& curves)
{
	const int numIter = 10;
	const double lambda = 0.5;
	for (int i = 0; i < numIter; ++i)
	{
		for (int j = 0; j < curves.size(); ++j)
		{
			MyPolyline newPoly = curves[j];
			for (int k = 1; k + 1 < curves[j].size(); ++k)
			{
				Eigen::Vector2d prev = curves[j][k - 1] - curves[j][k];
				Eigen::Vector2d next = curves[j][k + 1] - curves[j][k];
				double wPrev = 1 / prev.norm(), wNext = 1 / next.norm();
				double cosAngle = prev.normalized().dot(next.normalized());
				Eigen::Vector2d L = (wPrev*prev + wNext*next) / (wPrev + wNext);
				newPoly[k] += lambda*L;
			}
			curves[j] = newPoly;
		}
	}
}