#include "stdafx.h"
#include "Optimizer.h"
#include "TotalEnergy.h"
#include "LBFGS.h"
#include "Params.h"
#include <iostream>
#include <math.h>
#include <random>
#include "polynomial_energy.h"
#include "l2_regularizer.h"

//#define _USE_IPOPT_ 1

#ifdef _USE_IPOPT_
#include "IpIpoptApplication.hpp"
#include "MyNLP.hpp"
#endif


Eigen::VectorXcd optimize(cv::Mat & bwImg, const Eigen::MatrixXd & weight, const Eigen::MatrixXcd & tauNormalized, double beta, cv::Mat & mask, const Eigen::MatrixXi& indices)
{
	using namespace cv;
	int m = bwImg.rows;
	int n = bwImg.cols;
	int nnz = countNonZero(mask);
	std::cout << "nnz = " << nnz << std::endl;
	TotalEnergy fun(bwImg, weight, tauNormalized, beta, mask, indices,nnz);
	Eigen::VectorXd X(nnz * 4); //intial guess
	X.setZero();

#ifndef _USE_IPOPT_


	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-4;
	param.max_iterations = 2000;

	LBFGSpp::LBFGSSolver<double> solver(param);

	double fx;
	int niter = solver.minimize(fun, X, fx);

	std::cout << "Done in " << niter << " iterations" << std::endl;
	std::cout << "f(x) = " << fx << std::endl;

#else
	using namespace Ipopt;
	SmartPtr<TNLP> mynlp = new MyNLP(fun,X);

	// Create a new instance of IpoptApplication
	//  (use a SmartPtr, not raw)
	// We are using the factory, since this allows us to compile this
	// example with an Ipopt Windows DLL
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);

	// Change some options
	// Note: The following choices are only examples, they might not be
	//       suitable for your optimization problem.
	app->Options()->SetNumericValue("tol", 1e-2);
	app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	app->Options()->SetStringValue("linear_solver", "pardiso");
	app->Options()->SetNumericValue("print_level", 0);
	app->Options()->SetNumericValue("print_frequency_iter", 100);

	// The following overwrites the default name (ipopt.opt) of the
	// options file
	// app->Options()->SetStringValue("option_file_name", "hs071.opt");

	// Initialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) {
		std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
		return Eigen::VectorXcd();
	}

	// Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(mynlp);

	if (status == Solve_Succeeded) {
		std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
	}
	else {
		std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
	}
#endif
	Eigen::VectorXcd x_complex = X.head(X.size() / 2) + std::complex<double>(0, 1)*X.tail(X.size() / 2);
	return x_complex;
}

std::pair<Eigen::MatrixXd,Eigen::VectorXd> computeSmartWeights(const Eigen::MatrixXd& weight, cv::Mat& mask, const Eigen::MatrixXi& indices)
{
	int nnz = countNonZero(mask);
	int m = weight.rows();
	int n = weight.cols();
	Eigen::MatrixXd smartWeights(m, n);
	smartWeights.setOnes();
	for (int i = 1; i < m - 1; ++i)
		for (int j = 1; j < n - 1; ++j)
		{
			smartWeights(i, j) = 0.25 * ((1 - weight(i - 1, j)) + (1 - weight(i + 1, j)) + (1 - weight(i, j - 1)) + (1 - weight(i, j + 1)));
		}

	Eigen::VectorXd smartWeightsVector(nnz);

	for (int j = 0; j < n; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			if (mask.at<uchar>(i, j) == 0)
				continue;
			int idx = indices(i, j);
			smartWeightsVector(idx) = smartWeights(i, j);
		}
	}
	return std::make_pair(smartWeights, smartWeightsVector);
}

Eigen::VectorXcd optimizeByLinearSolve(cv::Mat& bwImg, const Eigen::MatrixXd& weight, const Eigen::MatrixXcd& tauNormalized, double beta, cv::Mat& mask, const Eigen::MatrixXi& indices)
{
	using namespace cv;
	int m = bwImg.rows;
	int n = bwImg.cols;
	int nnz = countNonZero(mask);
	std::cout << "nnz = " << nnz << std::endl;

	auto [A, b] = polynomial_energy_matrix(weight, tauNormalized, mask, indices);

	Eigen::MatrixXd onesMatrix;
	onesMatrix = Eigen::MatrixXd(m, n);
	onesMatrix.setOnes();
	Eigen::MatrixXcd g = tauNormalized * std::complex<double>(0, 1);
	auto [A2, b2] = polynomial_energy_matrix(onesMatrix, g, mask, indices);

	Eigen::SparseMatrix<double> L = laplacian_matrix(mask, indices, computeSmartWeights(weight,mask,indices).second);
	
	const double alpha = FRAME_FIELD_REGULARIZER_WEIGHT;
	Eigen::SparseMatrix<std::complex<double>> totalMatrix = 2 * (A + alpha*A2) + beta * 2.0 * L.cast<std::complex<double>>();
	Eigen::VectorXcd totalRhs = -2 * b.conjugate() - 2 * alpha * b2.conjugate();

	Eigen::ConjugateGradient<Eigen::SparseMatrix<std::complex<double>>, Eigen::Lower | Eigen::Upper> cg;
	cg.compute(totalMatrix);

	Eigen::VectorXcd result = cg.solve(totalRhs);
	if (cg.info() != Eigen::Success) {
		// solving failed
		return Eigen::VectorXcd();
	}
	else
		return result;
}

Eigen::VectorXcd optimizeByLinearSolve_holdingSomeFixed(cv::Mat& bwImg, const Eigen::MatrixXd& weight, const Eigen::MatrixXcd& tauNormalized, double beta, cv::Mat& mask, cv::Mat& fixedMask, const Eigen::MatrixXi& indices, const Eigen::MatrixXi& fixedIndices, const Eigen::VectorXcd& Xfixed)
{
	using namespace cv;
	int m = bwImg.rows;
	int n = bwImg.cols;
	int nnz = countNonZero(mask);
	int oldNnz = countNonZero(fixedMask);
	std::cout << "nnz = " << nnz << std::endl;

	Eigen::SparseMatrix<double> L = laplacian_matrix(mask, indices, computeSmartWeights(weight, mask, indices).second);
	//convert it to triplets
	std::vector<Eigen::Triplet<double>> newTrs;
	Eigen::VectorXcd rhs(nnz*2);
	rhs.setZero();

	std::map<int, int> newIdxToOldIdx;
	for (int j = 0; j < n; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			if ((mask.at<uchar>(i, j) == 0) || (fixedMask.at<uchar>(i,j)==0))
				continue;
			int idx = indices(i, j);
			int oldIdx = fixedIndices(i,j);
			newIdxToOldIdx.insert({ idx,oldIdx });
		}
	}

	auto oldIdx = [&newIdxToOldIdx](int idx)
	{
		auto it = newIdxToOldIdx.find(idx);
		if (it != newIdxToOldIdx.end())
			return it->second;
		else
			return -1;
	};

	std::vector<bool> fixEqnAdded(nnz); //says if we have already fixed that (new) index
	for (int k = 0; k < L.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
		{
			int newI = it.row(), newJ = it.col();
			if ((newI >= nnz) || (newJ >= nnz)) //L is (2n x 2n), for all variables
				continue;

			int oldI = oldIdx(newI), oldJ = oldIdx(newJ);
			if (oldI != -1)
			{
				if (!fixEqnAdded[newI])
				{
					fixEqnAdded[newI] = true;
					//fix x0
					newTrs.push_back({ newI, newI, 1.0 });
					rhs[newI] = Xfixed[oldI];

					//fix x2
					newTrs.push_back({ newI+nnz, newI+nnz, 1.0 });
					rhs[newI+nnz] = Xfixed[oldI+oldNnz];
				}
			}
			else if (oldJ != -1)
			{
				rhs(newI) -= it.value() * Xfixed(oldJ);
				rhs(newI+nnz) -= it.value() * Xfixed(oldJ+oldNnz);
			}
			else
			{
				newTrs.push_back({ (int)it.row(),(int)it.col(),it.value() }); //for x0
				newTrs.push_back({ (int)it.row()+nnz,(int)it.col()+nnz,it.value() }); //copy for x2
			}
		}



	//now let's fix the known part of the frame field
	

	Eigen::SparseMatrix<std::complex<double>> totalMatrix(nnz*2,nnz*2);
	totalMatrix.setFromTriplets(newTrs.begin(),newTrs.end());

	Eigen::ConjugateGradient<Eigen::SparseMatrix<std::complex<double>>, Eigen::Lower | Eigen::Upper> cg;
	cg.compute(totalMatrix);

	//std::ofstream fff("matdebug.m");
	//fff << "M = [" << Eigen::MatrixXd(totalMatrix.real()) << "];" << std::endl;
	//fff.close();

	Eigen::VectorXcd result = cg.solve(rhs);
	if (cg.info() != Eigen::Success) {
		// solving failed
		return Eigen::VectorXcd();
	}
	else
	{
		/*double error = 0;
		std::cout << "Fixed indices: ";
		for (auto it : newIdxToOldIdx)
		{
			std::complex<double> newVal, oldVal;
			std::cout << it.first << " ";
			newVal = result[it.first];
			oldVal = Xfixed[it.second];
			error += std::abs(std::pow(newVal - oldVal, 2));
		}*/
		return result;
	}
}
