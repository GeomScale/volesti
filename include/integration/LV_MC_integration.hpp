// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Integration and Optimization algorithm used here : https://www.cc.gatech.edu/~vempala/acg/www/www/papers/integration.pdf

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include "convex_bodies/hpolytope.h"
#include "Eigen/Eigen"
#include "generators/known_polytope_generators.h"
#include "boost_random_number_generator.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "ode_solvers/oracle_functors.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> HPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;

enum volumetype { CB ,CG ,SOB }; // Volume type for polytope


template
<
    typename EvaluationFunctor,
    typename GradientFunctor,
	typename Parameters,
    typename Polytope = HPOLYTOPE,
    typename Point,
    typename NT
>
NT lovasz_vempala_integrate(EvaluationFunctor &f,
                            GradientFunctor &grad_f,
							Parameters &params,
                            Polytope &P,
                            Point x0,
                            NT B,
                            unsigned int walk_length = 10,
                            NT epsilon = 0.1)
{
    unsigned int n = P.dimension();
    unsigned int m = (unsigned int) ceil(sqrt(n) * log(B));
    unsigned int k = (unsigned int) ceil(512 / pow(epsilon,2) * sqrt(n) * log(B));

    NT volume = volume_sequence_of_balls <BallWalk, RandomNumberGenerator, Polytope> (P, epsilon, walk_length);
    NT alpha_prev = (NT) 1 / B;
    NT alpha = (NT) 1 / B;
    NT log_W = log(volume);
    NT W_current = (NT) 0;

    RandomNumberGenerator rng(1);
	std::cerr << "n = " << n << " m = " << m << " k = " << k << " volume = " << volume << " log_W = " << log_W  << std::endl;
	std::cerr << "alpha = " << alpha << " alpha_prev = " << alpha_prev << " W_current = " << W_current << std::endl << std::endl;

    // Initialize HMC walks using EvaluationFunctor and GradientFunctor

	typedef LeapfrogODESolver<Point, NT, Polytope, GradientFunctor> Solver;

	HamiltonianMonteCarloWalk::parameters <NT, GradientFunctor> hmc_params(grad_f, n);

	HamiltonianMonteCarloWalk::Walk
	  <Point, Polytope, RandomNumberGenerator, GradientFunctor, EvaluationFunctor, Solver>
	  hmc(&P, x0, grad_f, f, hmc_params);

    // Check and evaluate for all samples breaks when variance > 1, i.e. alpha > 1
    for (int i = 1; i <= m && alpha < 1; i++ ) {

        alpha *= (1 + 1 / sqrt(n));
        W_current = 0;

        for (unsigned int j = 1; j <= k ; j++) {

            hmc.apply(rng, walk_length);
            W_current += exp(-f(hmc.x) * (alpha - alpha_prev));

        }

        W_current /= k;
        log_W += log(W_current);
        alpha_prev = alpha;
		params.update_temperature();
		std::cerr << "Value of alpha after i_th round " << alpha << std::endl;
    }

    return exp(log_W);    
}

/*
include/ode_solvers/oracle_functors 
To see how to define oracle_functors and gradient functors
https://github1s.com/GeomScale/volume_approximation/blob/develop/include/ode_solvers/oracle_functors.hpp

test/logconcave_sampling_test.cpp
https://github1s.com/GeomScale/volume_approximation/blob/develop/test/logconcave_sampling_test.cpp

HMC examples/logconcave
https://github.com/GeomScale/volume_approximation/tree/develop/examples/logconcave
*/

/*
  MT points = uniform sample points in ` n rows * k columns ` ( k is the number of points )
  RNG rng(1);    
  MT points(n,k); // n rows(dimension) * k columns(number of points)
  typename WalkType::template Walk <Polytope, RNG> walk(P, x0, rng);
      
  for (int i = 0; i < k; i++) {
      walk.apply(P, x0, walk_length, rng);
      points.col(i) = x0.getCoefficients();
  }
*/