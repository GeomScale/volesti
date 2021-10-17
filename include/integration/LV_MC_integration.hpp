// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Integration algorithm used here : https://www.cc.gatech.edu/~vempala/acg/www/www/papers/integration.pdf (Page-7)

#ifndef LOVASZ_VEMPALA_MC_INTEGRATION_HPP
#define LOVASZ_VEMPALA_MC_INTEGRATION_HPP

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

enum volumetype { CB,CG,SOB }; // Volume type for polytope

template
<
	typename EvaluationFunctor,
	typename GradientFunctor,
	typename WalkType,
	typename Polytope = HPOLYTOPE,
	typename Point,
	typename NT
>
NT lovasz_vempala_integrate(EvaluationFunctor &g,
							GradientFunctor &grad_g,
							Polytope &P,
							Point x0,
							NT beta = 1.0,
							volumetype voltype = SOB,
							unsigned int walk_length = 10,
							NT epsilon = 0.1)
{
	unsigned int n = P.dimension();
	NT B = 2 * n + 2 * log(1 / epsilon) + n * log(1 / beta);
	unsigned int m = (unsigned int) ceil(sqrt(n) * log(B));
	unsigned int k = 10000; // (unsigned int) ceil(512 / pow(epsilon,2) * sqrt(n) * log(B))/1000;

	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	OptimizationParameters opt_params(1, P.dimension(), g, grad_g);
	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	NT volume;

	switch (voltype) {
    case CB:     
        volume = volume_cooling_balls <AcceleratedBilliardWalk, RandomNumberGenerator, Polytope> (P, epsilon, walk_length).second; 
        break;
    case CG: 
        volume = volume_cooling_gaussians <GaussianBallWalk, RandomNumberGenerator, Polytope> (P, epsilon, walk_length);
        break;
    case SOB: 
        volume = volume_sequence_of_balls <AcceleratedBilliardWalk, RandomNumberGenerator, Polytope> (P, epsilon, walk_length);
        break;
    default:
        std::cerr << "Error in volume type: CB / SOB / CG" << std::endl;
        return -1;
    }

	std::cout << "Volume of the convex body = " << volume << std::endl;

	NT alpha = (NT) 1/B;
	NT alpha_prev = (NT) 0;
	NT log_W = log(volume);
	NT W_current = (NT) 0;

	RandomNumberGenerator rng(1);

	// Initialize HMC walks using EvaluationFunctor and GradientFunctor
	typedef LeapfrogODESolver<Point, NT, Polytope, NegativeGradientOptimizationFunctor> Solver;
	HamiltonianMonteCarloWalk::parameters <NT, NegativeGradientOptimizationFunctor> hmc_params(grad_f, n);
	HamiltonianMonteCarloWalk::Walk <Point, Polytope, RandomNumberGenerator, NegativeGradientOptimizationFunctor, NegativeLogprobOptimizationFunctor, Solver>
	  hmc(&P, x0, grad_f, f, hmc_params);

	// Burning samples for proper mixing
	typename WalkType::template Walk <Polytope, RandomNumberGenerator> walk(P, x0, rng);	  
	for (int i = 1; i <= k; i++) {
		walk.apply(P, x0, walk_length, rng);
	}

	// exit loop at alpha > 1
	bool loop_run = true;
	while (loop_run) {
		
		if (alpha > 1) { alpha = 1; loop_run = false; }
		opt_params.set_temperature(alpha_prev);
		W_current = 0;

		for (unsigned int j = 1; j <= k ; j++) {

			hmc.apply(rng, walk_length);
			W_current += exp(-g(hmc.x) * (alpha - alpha_prev));
			
		}
		
		W_current /= k;
		log_W += log(W_current);
		alpha_prev = alpha;
		alpha *= (1 + 1 / sqrt(n));		
	}

	return exp(log_W);    
}

#endif
