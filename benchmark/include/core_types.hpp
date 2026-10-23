#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <boost/random.hpp>
#include "Eigen/Eigen"

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp" 

// In this file you can find the core type definitions and commonly used aliases
// for the project.
// 
// It defines:
// Numeric type (NT) used across computations.
// Geometric kernel (Cartesian) and Point representation.
// Matrix (MT) and vector (VT) types using Eigen (similar to NumPy in C++).
// Random number generator type based on Boost.
// HPolytope type used to represent convex bodies.

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef HPolytope<Point> HPOLYTOPE;

//it is only initialized once inside src/benchmark_utils.cpp
extern PushBackWalkPolicy push_back_policy;