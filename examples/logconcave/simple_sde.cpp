// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include <string>
#include <typeinfo>
#include <chrono>

#include "Eigen/Eigen"

#include "ode_solvers.hpp"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/known_polytope_generators.h"

template <typename NT>
void run_main(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
    typedef IsotropicQuadraticFunctor::GradientFunctor<Point> func;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;

    unsigned int dim = 1;

    // Random number generator
    RandomNumberGenerator rng(1);

    // Define functor parameters (order, multplier etc.)
    IsotropicQuadraticFunctor::parameters<NT> params;
    params.order = 2;
    params.alpha = NT(1);

    // Create functor
	func F(params);

    // Initial conditions
    Point x0 = Point(dim);
    Point v0 = Point::all_ones(dim);
    pts q{x0, v0};

    // Create solver object
    RandomizedMipointSDESolver<Point, NT, Hpolytope, func, RandomNumberGenerator>
    randomized_midpoint_solver =
      RandomizedMipointSDESolver<Point, NT, Hpolytope, func, RandomNumberGenerator>
      (0, 0.01, q, F, bounds{NULL, NULL}, 1.0 / F.params.L);

    // Simulate the ODE for 1000 steps and log the state (x, v) to stdout
    int n_steps = 1000;

    for (int i = 0; i < n_steps; i++) {
      randomized_midpoint_solver.step(rng);
      randomized_midpoint_solver.print_state();
    }

}

int main(void) {
  run_main<double>();
  return 0;
}
