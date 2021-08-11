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
#include "LV_MC_optimization.hpp"
#include "ode_solvers/oracle_functors.hpp"
#include "generators/known_polytope_generators.h"
#include "boost_random_number_generator.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "misc.h"

typedef double NT;

template <typename NT>
void call_test_lv_optimization() {
	
	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
	typedef typename HPOLYTOPE::MT MT;
    typedef typename HPOLYTOPE::VT VT;

	typedef GaussianFunctor::FunctionFunctor <MT,NT,Point> EvaluationFunctor;
	typedef GaussianFunctor::GradientFunctor <MT,NT,Point> GradientFunctor;
    typedef GaussianFunctor::parameters <MT,NT,Point> Parameters;

    HPOLYTOPE HP = generate_cube <HPOLYTOPE> (3, false);
    unsigned int n = HP.dimension();
	std::pair <Point,NT> inner_ball = HP.ComputeInnerBall();
	Point x0 = inner_ball.first;

    Parameters gaussian_params(x0, 1.0);
    EvaluationFunctor g(gaussian_params);
    GradientFunctor grad_g(gaussian_params);

    NT beta = 1.0;

    std::pair<Point,NT> pair_functor= lovasz_vempala_optimize
      <EvaluationFunctor, GradientFunctor, Parameters, BilliardWalk, HPOLYTOPE, Point, MT, NT>
      (g, grad_g, gaussian_params, HP, x0);

}

TEST_CASE("optimize") {
    call_test_lv_optimization<double>();
}