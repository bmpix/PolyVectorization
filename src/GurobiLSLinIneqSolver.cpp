#include "stdafx.h"
#include "GurobiLSLinIneqSolver.h"
#include <gurobi_c++.h>

class GurobiLSLinIneqSolver::Impl
{
	// Nothing yet!
};

GurobiLSLinIneqSolver::GurobiLSLinIneqSolver()
	: pImpl_(new Impl)
{
	// Nothing yet!
}

GurobiLSLinIneqSolver::~GurobiLSLinIneqSolver()
{
	// Nothing yet!
}

template<typename T>
void IterateSparseMatrixNNZ(const Eigen::SparseMatrix<double>& mat, T func)
{
	for (int k=0; k<mat.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
		{
			//it.value();
			//it.row();   // row index
			//it.col();   // col index (here it is equal to k)
			//it.index(); // inner index, here it is equal to it.row()
			func(it.row(), it.col(), it.value());
		}
	}
}

std::vector<double> GurobiLSLinIneqSolver::solve( const Eigen::SparseMatrix<double>& H, const std::vector<double>& f, const Eigen::SparseMatrix<double>& A, const std::vector<double> b, const Eigen::SparseMatrix<double>& Aeq, const std::vector<double>& bEq, int nUnknowns )
{
	return solve(H,f,A,b,Aeq,bEq,nUnknowns,-1.0,1.0);
}

std::vector<double> GurobiLSLinIneqSolver::solve( const Eigen::SparseMatrix<double>& H, const std::vector<double>& f, const Eigen::SparseMatrix<double>& A, const std::vector<double> b, const Eigen::SparseMatrix<double>& Aeq, const std::vector<double>& bEq, int nUnknowns, double lb, double ub )
{
	try
	{
		auto env = std::unique_ptr<GRBEnv>(new GRBEnv());
		GRBModel model(*env);
		/*if (Inbetweening::Settings::getInstance().debugLevel() <= 2)
		{
			env->set(GRB_IntParam_OutputFlag, 0);
		}*/

		// Define variable properties
		auto lbs = std::unique_ptr<double[]>(new double[nUnknowns]);
		auto ubs = std::unique_ptr<double[]>(new double[nUnknowns]);
		auto varTypes = std::unique_ptr<char[]>(new char[nUnknowns]);
		for (int i = 0; i < nUnknowns; ++i)
		{
			lbs[i] = lb;
			ubs[i] = ub;
			varTypes[i] = GRB_CONTINUOUS;
		}

		// Add variables to the model
		auto vars = std::unique_ptr<GRBVar[]>(model.addVars(lbs.get(), ubs.get(), NULL, varTypes.get(), NULL, nUnknowns));
		model.update();

		// Populate objective
		{
			GRBQuadExpr obj = 0;

			// Quadratic component
			assert(nUnknowns == H.rows() && nUnknowns == H.cols());
			IterateSparseMatrixNNZ(H, [&obj,&vars] (int row, int col, double val) {
				obj += (0.5 * val) * vars[row] * vars[col];
			});

			// Linear component
			assert(nUnknowns == f.size());
			for (int i = 0; i < f.size(); ++i)
			{
				if (f[i] != 0.0)
					obj += f[i] * vars[i];
			}

			model.setObjective(obj, GRB_MINIMIZE);
		}

		// Populate Equality constraints
		if (Aeq.rows() > 0)
		{
			const int numeq = Aeq.rows();
			assert(numeq == bEq.size());
			//std::cout << "Gurobi populating " << numeq << " equality constraints." << std::endl;

			auto eq = std::unique_ptr<GRBLinExpr[]>(new GRBLinExpr[numeq]);
			IterateSparseMatrixNNZ(Aeq, [&eq,&vars] (int row, int col, double val) {
				eq[row] += val * vars[col];
			});
					
			for (int i=0; i < numeq; ++i)
				model.addConstr(eq[i], GRB_EQUAL, bEq[i]);
		}

		// Populate Inequality constraints
		if (A.rows() > 0)
		{
			const int numineq = A.rows();
			assert(numineq == b.size());
			//std::cout << "Gurobi populating " << numineq << " inequality constraints." << std::endl;

			auto ineq = std::unique_ptr<GRBLinExpr[]>(new GRBLinExpr[numineq]);

			IterateSparseMatrixNNZ(A, [&ineq,&vars] (int row, int col, double val) {
				ineq[row] += val * vars[col];
			});

			for (int i=0; i < numineq; ++i)
				model.addConstr(ineq[i], GRB_LESS_EQUAL, b[i]);
		}

		model.update();
		/*if (Inbetweening::Settings::getInstance().debugLevel() > 3)
		{
			model.write("GurobiSolve.lp");
		}*/

		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL
			|| model.get(GRB_IntAttr_Status) == GRB_SUBOPTIMAL)
		{
			double objvalP = model.get(GRB_DoubleAttr_ObjVal);

			//std::cout << "Gurobi found optimal solution: objective function value " << objvalP << std::endl;

			std::vector<double> solution(nUnknowns, 0.0);
			for (int i = 0; i < nUnknowns; i++)
				solution[i] = vars[i].get(GRB_DoubleAttr_X);
		
			return solution;
		}
	}
	catch(GRBException e)
	{
		std::cout << "Gurobi Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch(...)
	{
		std::cout << "Exception during Gurobi optimization!" << std::endl;
	}

	return std::vector<double>();
}

std::vector<double> GurobiLSLinIneqSolver::solveLeastSquares(
	const std::vector<Eigen::Triplet<double>>& A, const std::vector<double>& b,
	const std::vector<Eigen::Triplet<double>>& A_inequality, const std::vector<double>& b_inequality,
	const std::vector<Eigen::Triplet<double>>& A_equality, const std::vector<double>& b_equality,
	int nUnknowns, double lb, double ub )
{
	try
	{
		auto env = std::unique_ptr<GRBEnv>(new GRBEnv());
		/*if (Inbetweening::Settings::getInstance().debugLevel() <= 2)
		{
			env->set(GRB_IntParam_OutputFlag, 0);
		}*/
		env->set(GRB_DoubleParam_BarConvTol, 1e-6);

		GRBModel model(*env);

		// Define variable properties
		auto lbs = std::unique_ptr<double[]>(new double[nUnknowns]);
		auto ubs = std::unique_ptr<double[]>(new double[nUnknowns]);
		auto varTypes = std::unique_ptr<char[]>(new char[nUnknowns]);
		for (int i = 0; i < nUnknowns; ++i)
		{
			lbs[i] = lb;
			ubs[i] = ub;
			varTypes[i] = GRB_CONTINUOUS;
		}

		// Add variables to the model
		auto vars = std::unique_ptr<GRBVar[]>(model.addVars(lbs.get(), ubs.get(), NULL, varTypes.get(), NULL, nUnknowns));
		model.update();

		// Populate objective
		{
			const int numRows = b.size();

			// Individual errors
			auto linexp = std::unique_ptr<GRBLinExpr[]>(new GRBLinExpr[numRows]);
			for (auto& tri : A)
			{
				assert(tri.row() < numRows);
				assert(tri.col() < nUnknowns);
				linexp[tri.row()] += tri.value() * vars[tri.col()];
			}

			// Construct total quadratic expression
			GRBQuadExpr quadexp = 0;
			for (int i = 0; i < numRows; ++i)
			{
				linexp[i] -= b[i];
				quadexp += linexp[i] * linexp[i];
			}

			model.setObjective(quadexp, GRB_MINIMIZE);
		}

		// Populate Equality constraints
		if (b_equality.size() > 0)
		{
			const int numeq = b_equality.size();
			//std::cout << "Gurobi populating " << numeq << " equality constraints." << std::endl;
			
			auto eq = std::unique_ptr<GRBLinExpr[]>(new GRBLinExpr[numeq]);
			for (auto& tri : A_equality)
			{
				assert(tri.row() < numeq);
				assert(tri.col() < nUnknowns);
				eq[tri.row()] += tri.value() * vars[tri.col()];
			}
			
			for (int i=0; i < numeq; ++i)
				model.addConstr(eq[i], GRB_EQUAL, b_equality[i]);
		}

		// Populate Inequality constraints
		if (b_inequality.size() > 0)
		{
			const int numineq = b_inequality.size();
			//std::cout << "Gurobi populating " << numineq << " inequality constraints." << std::endl;

			auto ineq = std::unique_ptr<GRBLinExpr[]>(new GRBLinExpr[numineq]);
			for (auto& tri : A_inequality)
			{
				assert(tri.row() < numineq);
				assert(tri.col() < nUnknowns);
				ineq[tri.row()] += tri.value() * vars[tri.col()];
			}

			for (int i=0; i < numineq; ++i)
				model.addConstr(ineq[i], GRB_LESS_EQUAL, b_inequality[i]);
		}

		model.update();
		/*if (Inbetweening::Settings::getInstance().debugLevel() > 3)
		{
			model.write("GurobiSolve_LeastSquares.lp");
		}*/

		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL
			|| model.get(GRB_IntAttr_Status) == GRB_SUBOPTIMAL)
		{
			double objvalP = model.get(GRB_DoubleAttr_ObjVal);

			//std::cout << "Gurobi found optimal solution: objective function value " << objvalP << std::endl;

			std::vector<double> solution(nUnknowns, 0.0);
			for (int i = 0; i < nUnknowns; i++)
				solution[i] = vars[i].get(GRB_DoubleAttr_X);

			return solution;
		}
	}
	catch(GRBException e)
	{
		std::cout << "Gurobi Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch(...)
	{
		std::cout << "Exception during Gurobi optimization!" << std::endl;
	}

	return std::vector<double>();
}
