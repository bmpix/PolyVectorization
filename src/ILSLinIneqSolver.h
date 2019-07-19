#ifndef _LS_LIN_INEQ_SOLVE_H_
#define _LS_LIN_INEQ_SOLVE_H_
#include <Eigen/Core>
#include <Eigen/Sparse>

/**
 * This is an interface header for solving Linear (Overdetermined) System with linear inequality constraints.
 * Minimizes 1/2*x'*H*x + f'*x subject to the restrictions A*x ≤ b, additional restrictions Aeq*x = beq.
 */
class ILSLinIneqSolver
{
public:
	typedef Eigen::Triplet<double> Triplet;

	/**
	 * Main solving interface.
	 */
	virtual std::vector<double> solve(const Eigen::SparseMatrix<double>& H, const std::vector<double>& f, const Eigen::SparseMatrix<double>& A, const std::vector<double> b, const Eigen::SparseMatrix<double>& Aeq, const std::vector<double>& bEq, int nUnknowns) = 0;
};

#endif // _LS_LIN_INEQ_SOLVE_H_
