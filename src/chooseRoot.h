#ifndef _CHOOSE_ROOTS_H_
#define _CHOOSE_ROOTS_H_

#include "typedefs.h"

Eigen::Vector2d chooseRoot(std::complex<double> root1, std::complex<double> root2, Eigen::Vector2d dir, Eigen::Vector2d* otherRoot = nullptr);

#endif
