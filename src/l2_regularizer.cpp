#include "stdafx.h"
#include "l2_regularizer.h"
#include <map>
#include <array>
#include <numeric>
#include "typedefs.h"
//#define doesNeighborExist(i,j,m,n,vertical,left)  ((left)? ((vertical)?((i)>0):((j)>0)) : (((vertical)?(i):(j))<((vertical)?(m):(n))-1) )

std::pair<double, Eigen::VectorXcd> l2_regularizer(const Eigen::VectorXcd & X, const Eigen::MatrixXd & D, const cv::Mat & mask, Eigen::MatrixXi& indices, bool computeGrad)
{
	int nonzeros = X.size()/2;
	assert(nonzeros == countNonZero(mask));

	double energy = 0;
	Eigen::VectorXcd grad;
	if (computeGrad)
	{
		grad = Eigen::VectorXcd(nonzeros * 2);
		grad.setZero();
	}

	Eigen::VectorXcd dx(nonzeros);
	
	int m = mask.rows, n = mask.cols;
	for (int k = 0; k < 2; ++k) //x0 or x2
	{
		for (int dir = 0; dir < 2; dir++) //horizontal or vertical
		{
			std::vector<double> energies(nonzeros,0);
		#pragma omp parallel for
			for (int j = 0; j < n; ++j)
			{
				for (int i = 0; i < m; ++i)
				{
					if (mask.at<uchar>(i, j) == 0)
						continue;
					int idx = indices(i, j);
					std::pair<int, int> leftNeighbor, rightNeighbor;
					bool useLeft = useNeighbor(i, j, m, n, (dir == 1), true, mask, leftNeighbor);
					bool useRight = useNeighbor(i, j, m, n, (dir == 1), false, mask, rightNeighbor);

					if (useLeft & useRight)
						dx(idx) = 0.5*(X(indices(rightNeighbor.first, rightNeighbor.second) + k*nonzeros) - X(indices(leftNeighbor.first, leftNeighbor.second) + k*nonzeros));
					else if (useLeft)
						dx(idx) = X(idx + k*nonzeros) - X(indices(leftNeighbor.first, leftNeighbor.second) + k*nonzeros);
					else if (useRight)
						dx(idx) = X(indices(rightNeighbor.first, rightNeighbor.second) + k*nonzeros) - X(idx + k*nonzeros);
					else
						dx(idx) = 0;

					auto dx1 = dx(idx);

					//energy += D(i, j) * (dx1.real()*dx1.real() + dx1.imag()*dx1.imag());
					energies[idx] += D(i, j) * (dx1.real()*dx1.real() + dx1.imag()*dx1.imag());
				}
			}

			double curEnergy = std::accumulate(energies.begin(), energies.end(), 0.0);
			energy += curEnergy;

			if (computeGrad)
			{
#pragma omp parallel for
				for (int i = 0; i < m; ++i)
				{
					for (int j = 0; j < n; ++j)
					{
						if (mask.at<uchar>(i, j) == 0)
							continue;

						int idx = indices(i, j);
						for (int sign : {-1, 1}) //leftRight
						{
							std::pair<int, int> neighbor;
							bool shouldIUseNeighbor = useNeighbor(i, j, m, n, (dir == 1), sign == 1, mask, neighbor);
							bool shouldIUseNeighborsNeighbor = false;
							std::pair<int, int> tmpNeighbor;
							if (doesNeighborExist(i, j, m, n, (dir == 1), sign == 1))
								shouldIUseNeighborsNeighbor = useNeighbor(neighbor.first, neighbor.second, m, n, (dir == 1), sign == 1, mask, tmpNeighbor);

							if (shouldIUseNeighbor)
							{
								int myNeighbor = indices(neighbor.first, neighbor.second);
								if (shouldIUseNeighborsNeighbor)
									grad(idx + k*nonzeros) += double(sign)*dx(myNeighbor)*D(neighbor.first, neighbor.second);
								else
									grad(idx + k*nonzeros) += 2 * double(sign)*dx(myNeighbor)*D(neighbor.first, neighbor.second);
							}
							else
							{
								grad(idx + k*nonzeros) -= 2 * double(sign)*dx(idx)*D(i, j);
							}
						}

					}
				}
			}
		}
	}

	return std::make_pair(energy, grad);
}
