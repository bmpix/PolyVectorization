#ifndef _EQUATIONS_H_
#define _EQUATIONS_H_
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> Triplet;
struct Equations
{
	std::vector<Triplet> ineq, eq, exactEq;
	std::vector<double> ineqB, eqB, exactB;

	void add(Equations other)
	{
		auto offsetRow = [](std::vector<Triplet>& tripletsInGlobalIndices, int offset)
		{
			for (auto& tr : tripletsInGlobalIndices)
				tr = Triplet(tr.row() + offset, tr.col(), tr.value());
		};

		offsetRow(other.ineq, ineqB.size());
		ineq.insert(ineq.end(), other.ineq.begin(), other.ineq.end());
		ineqB.insert(ineqB.end(), other.ineqB.begin(), other.ineqB.end());

		offsetRow(other.eq, eqB.size());
		eq.insert(eq.end(), other.eq.begin(), other.eq.end());
		eqB.insert(eqB.end(), other.eqB.begin(), other.eqB.end());

		offsetRow(other.exactEq, exactB.size());
		exactEq.insert(exactEq.end(), other.exactEq.begin(), other.exactEq.end());
		exactB.insert(exactB.end(), other.exactB.begin(), other.exactB.end());
	}

	void addEqLHS(int index, double value)
	{
		eq.push_back(Triplet(eqB.size(), index, value));
	}

	void addExactEqLHS(int index, double value)
	{
		exactEq.push_back(Triplet(exactB.size(), index, value));
	}

	void addEqRHSAndFinishEquation(double rhs)
	{
		eqB.push_back(rhs);
	}

	void addExactRHSAndFinishEquation(double rhs)
	{
		exactB.push_back(rhs);
	}

	void addIneqLHS(int index, double value)
	{
		ineq.push_back(Triplet(ineqB.size(), index, value));
	}

	void addIneqRHSAndFinishEquation(double rhs)
	{
		ineqB.push_back(rhs);
	}
};

#endif