// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Integration and Optimization algorithm used here : https://www.cc.gatech.edu/~vempala/acg/www/www/papers/integration.pdf

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
typedef std::pair<Point, NT> PairFunctor;

enum volumetype { CB,CG,SOB }; // Volume type for polytope

template
<
	typename EvaluationFunctor,
	typename GradientFunctor,
	typename Parameters,
	typename WalkType,
	typename Polytope = HPOLYTOPE,
	typename Point,
	typename NT
>
NT lovasz_vempala_integrate(EvaluationFunctor &g,
							GradientFunctor &grad_g,
							Parameters &params,
							Polytope &P,
							Point x0,
							NT B,
							volumetype voltype = SOB,
							unsigned int walk_length = 10,
							NT epsilon = 0.1)
{
	unsigned int n = P.dimension();
	unsigned int m = (unsigned int) ceil(sqrt(n) * log(B));
	unsigned int k = (unsigned int) ceil(512 / pow(epsilon,2) * sqrt(n) * log(B));

	NT volume = 0;

	switch (voltype) {
    case CB:     
        volume = volume_cooling_balls <BallWalk, RandomNumberGenerator, Polytope> (P, epsilon, walk_length).second; 
        break;
    case CG: 
        volume = volume_cooling_gaussians <GaussianBallWalk, RandomNumberGenerator, Polytope> (P, epsilon, walk_length);
        break;
    case SOB: 
        volume = volume_sequence_of_balls <BallWalk, RandomNumberGenerator, Polytope> (P, epsilon, walk_length);
        break;
    default:
        std::cerr << "Error in volume type: CB / SOB / CG" << std::endl;
        return -1;
    }

	NT alpha_prev = (NT) 1 / B;
	NT alpha = (NT) 1 / B;
	NT log_W = log(volume);
	NT W_current = (NT) 0;

	RandomNumberGenerator rng(1);

	// Burning samples for proper mixing
	typename WalkType::template Walk <Polytope, RandomNumberGenerator> walk(P, x0, rng);	  
	for (int i = 1; i <= k; i++) {
		walk.apply(P, x0, walk_length, rng);
	}
	std::cout << "Print x0: " ; x0.print();

	std::cerr << "B = " << B << " n = " << n << " m = " << m << " k = " << k << " volume = " << volume << " log_W = " << log_W  << std::endl;
	std::cerr << "alpha = " << alpha << " alpha_prev = " << alpha_prev << " W_current = " << W_current << std::endl << std::endl;

	// Initialize HMC walks using EvaluationFunctor and GradientFunctor

	typedef LeapfrogODESolver<Point, NT, Polytope, GradientFunctor> Solver;

	HamiltonianMonteCarloWalk::parameters <NT, GradientFunctor> hmc_params(grad_g, n);

	HamiltonianMonteCarloWalk::Walk
		<Point, Polytope, RandomNumberGenerator, GradientFunctor, EvaluationFunctor, Solver>
		hmc(&P, x0, grad_g, g, hmc_params);

	// Check and evaluate for all samples breaks when variance > 1, i.e. alpha > 1
	while(alpha <= 1) {

		alpha *= (1 + 1 / sqrt(n));
		params.set_temperature(alpha);
		W_current = 0;

		for (unsigned int j = 1; j <= k ; j++) {

			hmc.apply(rng, walk_length);
			W_current += exp(-g(hmc.x) * (alpha - alpha_prev));
			// std::cout << hmc.x.getCoefficients().transpose() << std::endl;
			
		}

		W_current /= k;
		log_W += log(W_current);
		alpha_prev = alpha;
		std::cerr << "After i_th round | alpha = " << alpha << " | exp(log_W) = " << exp(log_W) << std::endl;
	}

	return exp(log_W);    
}

template
<
	typename EvaluationFunctor,
	typename GradientFunctor,
	typename Parameters,
	typename WalkType,
	typename Polytope = HPOLYTOPE,
	typename Point,
	typename NT
>
std::pair<Point, NT> lovasz_vempala_optimize(EvaluationFunctor &g,
											 GradientFunctor &grad_g,
											 Parameters &params,
											 Polytope &P,
											 Point x0,
											 NT B,
											 NT C0 = 1,
											 unsigned int walk_length = 10,
											 NT delta = 1,
											 NT epsilon = 0.1)
{
	unsigned int n = P.dimension();
	unsigned int m = (unsigned int)	ceil(sqrt(n) * log(2 * B * (n + log(1 / delta)) / epsilon ));
	unsigned int k = (unsigned int) ceil(C0 * n * pow(log(n), 5));

	RandomNumberGenerator rng(1);
	
	// Burning samples for proper mixing
	typename WalkType::template Walk <Polytope, RandomNumberGenerator> walk(P, x0, rng);	  
	for (int i = 1; i <= k; i++) {
		walk.apply(P, x0, walk_length, rng);
	}

	// Initialize HMC walks using EvaluationFunctor and GradientFunctor

	typedef LeapfrogODESolver<Point, NT, Polytope, GradientFunctor> Solver;

	HamiltonianMonteCarloWalk::parameters <NT, GradientFunctor> hmc_params(grad_g, n);

	HamiltonianMonteCarloWalk::Walk
		<Point, Polytope, RandomNumberGenerator, GradientFunctor, EvaluationFunctor, Solver>
		hmc(&P, x0, grad_g, g, hmc_params);

	
	NT alpha = 1 / B;
	Point calculated_point = x0, max_point = x0;
	NT calculated_value = exp(-g(hmc.x)), max_value = calculated_value;

	// Check and evaluate for all samples breaks when variance > 1, i.e. alpha > 1
	for (int i = 1; i <= m ; i++) {

		alpha *= (1 + 1 / sqrt(n));
		params.set_temperature(alpha);

		for (unsigned int j = 1; j <= k ; j++) {

			hmc.apply(rng, walk_length);
			calculated_value = exp(-g(hmc.x));
			if (calculated_value > max_value) {
				max_point = hmc.x;
				max_value = calculated_value;
			}
			
		}
	
	}

	std::pair<Point, NT> values(max_point, max_value);
	// values.insert(std::pair<int, int>(0, 42));    
	std::cout << "Max point check = " ; values.first.print();
	std::cout << "Max value check = " << values.second << std::endl;

	return values;  // (Point)values.first = X s.t. max_f = max(exp(-g(X))), (NT)values.second = max_f

}

#endif

/*
include/ode_solvers/oracle_functors 
To see how to define oracle_functors and gradient functors
https://github1s.com/GeomScale/volume_approximation/blob/develop/include/ode_solvers/oracle_functors.hpp

test/logconcave_sampling_test.cpp
https://github1s.com/GeomScale/volume_approximation/blob/develop/test/logconcave_sampling_test.cpp

HMC examples/logconcave
https://github.com/GeomScale/volume_approximation/tree/develop/examples/logconcave
*/
