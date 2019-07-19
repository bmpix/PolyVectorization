#ifndef _GUROBI_LS_LIN_INEQ_SOLVE_H_
#define _GUROBI_LS_LIN_INEQ_SOLVE_H_

#include "ILSLinIneqSolver.h"
#include <memory>

/**
 * This is an interface header for solving Linear (Overdetermined) System with linear inequality constraints.
 * Minimizes 1/2*x'*H*x + f'*x subject to the restrictions A*x ≤ b, additional restrictions Aeq*x = beq.
 */
class GurobiLSLinIneqSolver : public ILSLinIneqSolver
{
public:
	GurobiLSLinIneqSolver();
	~GurobiLSLinIneqSolver();

	/**
	 * Main solving interface.
	 * Assumes that lower bound is -1.0, upper bound is 1.0.
	 */
	virtual std::vector<double> solve(const Eigen::SparseMatrix<double>& H, const std::vector<double>& f, const Eigen::SparseMatrix<double>& A, const std::vector<double> b, const Eigen::SparseMatrix<double>& Aeq, const std::vector<double>& bEq, int nUnknowns);

	/**
	 * Main solving interface with custom lower and upper bounds.
	 */
	std::vector<double> solve(const Eigen::SparseMatrix<double>& H, const std::vector<double>& f, const Eigen::SparseMatrix<double>& A, const std::vector<double> b, const Eigen::SparseMatrix<double>& Aeq, const std::vector<double>& bEq, int nUnknowns, double lb, double ub);

	/**
	 * Solving interface using triplets. Below is a simple example.
	 * Note that Gurobi requires you to specify a lower and upper bound for each variable.
	 * If you don't specify it, Gurobi will assume the lower bound is 0.0.
	 *
	 * 	// Test simple solve
	 * 	GurobiLSLinIneqSolver temp;
	 * 	auto testTriplets = emptyTriplets;
	 * 	auto testVec = emptyVec;
	 * 
	 * 	testTriplets.push_back(Triplet(0, 0, 1.0));
	 * 	testTriplets.push_back(Triplet(1, 0, 1.0));
	 * 	testTriplets.push_back(Triplet(2, 1, 1.0));
	 * 	testTriplets.push_back(Triplet(3, 1, 1.0));
	 * 	testVec.push_back(1.0);
	 * 	testVec.push_back(1.0);
	 * 	testVec.push_back(2.0);
	 * 	testVec.push_back(3.0);
	 * 
	 * 	auto X_gurobi_LS = temp.solveLeastSquares(testTriplets, testVec, emptyTriplets, emptyVec, emptyTriplets, emptyVec, 2, 0.0, 100.0);
	 * 	std::cout << "Gurobi test solve: " << std::endl;
	 * 	for (auto& entry: X_gurobi_LS)
	 * 		std::cout << entry << " ";
	 * 	std::cout << std::endl;
	 *
	 * --> Should output 1.0 and 2.5
	 */
	std::vector<double> solveLeastSquares( const std::vector<Eigen::Triplet<double>>& A, const std::vector<double>& b, const std::vector<Eigen::Triplet<double>>& A_inequality, const std::vector<double>& b_inequality, const std::vector<Eigen::Triplet<double>>& A_equality, const std::vector<double>& b_equality, int nUnknowns, double lb, double ub );

private:
	// Disable copying of this class
	GurobiLSLinIneqSolver(const GurobiLSLinIneqSolver& t);
	GurobiLSLinIneqSolver& operator=(const GurobiLSLinIneqSolver& t);

	class Impl;
	std::unique_ptr<Impl> pImpl_;
};

#endif // _GUROBI_LS_LIN_INEQ_SOLVE_H_
