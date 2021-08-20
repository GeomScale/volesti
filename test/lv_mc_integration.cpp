// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "convex_bodies/hpolytope.h"
#include "Eigen/Eigen"
#include "doctest.h"
#include "LV_MC_integration.hpp"
#include "ode_solvers/oracle_functors.hpp"
#include "generators/known_polytope_generators.h"
#include "boost_random_number_generator.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "misc.h"

typedef double NT;

template <typename NT>
void test_values (NT computed, NT expected, NT exact) {
	std::cout << "Computed integration value = " << computed << std::endl;
	std::cout << "Expected integration value = " << expected << std::endl;
	std::cout << "Exact integration value = " << exact << std::endl;
	std::cout << "Relative error (expected) = " << std::abs((computed - expected)/expected) << std::endl;
	std::cout << "Relative error (exact) = " << std::abs((computed - exact)/exact) << std::endl ;
	CHECK(((std::abs((computed - expected)/expected) < 0.00001) || (std::abs((computed - exact)/exact) < 0.2)));
}

template <typename NT>
void cubes_test_lovasz_vempala_integrate (const int dimension, NT expected, NT exact, NT beta = 1.0, volumetype voltype = SOB) { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef IsotropicQuadraticFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef IsotropicQuadraticFunctor::GradientFunctor <Point> GradientFunctor;
	typedef IsotropicQuadraticFunctor::parameters<NT> IsoParameters;

	IsoParameters params;
	params.alpha = (NT) 2;
	GradientFunctor grad_g(params);
	EvaluationFunctor g(params);

	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (dimension, false);
	unsigned int n = HP.dimension();;
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();;
	Point x0 = inner_ball.first;

	NT integral_value = lovasz_vempala_integrate 
	  <EvaluationFunctor, GradientFunctor, IsoParameters, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
	  (g, grad_g, params, HP, x0, beta, voltype, 5, 0.1);
	
	test_values(integral_value, expected, exact);
}

template<typename NT>
void call_cubes_test_lovasz_vempala_integrate() {
	cubes_test_lovasz_vempala_integrate <NT> (1, 1.50, 1.49364, 1.0, CB);
	cubes_test_lovasz_vempala_integrate <NT> (2, 2.24, 2.23);
	cubes_test_lovasz_vempala_integrate <NT> (3, 3.34, 3.33225);
	cubes_test_lovasz_vempala_integrate <NT> (5, 7.50, 7.4325);
	cubes_test_lovasz_vempala_integrate <NT> (7, 16.77, 16.5852);
}

TEST_CASE("cubes") {
	call_cubes_test_lovasz_vempala_integrate<double>();
}
