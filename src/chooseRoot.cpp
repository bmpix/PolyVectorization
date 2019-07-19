#include "stdafx.h"
#include "chooseRoot.h"
#include <iostream>

Eigen::Vector2d chooseRoot(std::complex<double> root1, std::complex<double> root2, Eigen::Vector2d dir, Eigen::Vector2d * otherRoot)
{
	Eigen::Vector2d eigRoot1 = _toEig(root1), eigRoot2 = _toEig(root2);
	eigRoot1.normalize(); eigRoot2.normalize();
	std::array<Eigen::Vector2d, 4> roots = { eigRoot1, -eigRoot1, eigRoot2, -eigRoot2 };
	
	std::array<double, 4> dots;
	for (int i = 0; i < 4; ++i)
		dots[i] = roots[i].dot(dir);

	auto it = std::max_element(dots.begin(), dots.end());
	int idx = it - dots.begin();

	if (otherRoot)
		*otherRoot = roots[(idx + 2) % 4];

	if (fabs(dots[idx] - dots[(idx + 2) % 4]) < 1e-3)
		return 0.5*(roots[idx] + roots[(idx + 2) % 4]);

	return roots[idx]; 
}
