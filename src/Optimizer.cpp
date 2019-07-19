#include "stdafx.h"
#include "Optimizer.h"
#include "TotalEnergy.h"
#include "LBFGS.h"
#include <iostream>
#include <math.h>
#include <random>

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
