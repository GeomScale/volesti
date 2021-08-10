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

enum volumetype { CB,CG,SOB }; // Volume type for polytope

template
<
	typename EvaluationFunctor,
	typename GradientFunctor,
	typename Parameters,
	typename WalkType,
	typename Polytope,
	typename Point,
	typename NT
>
NT lovasz_vempala_integrate(EvaluationFunctor &g,
							GradientFunctor &grad_g,
							Parameters &params,
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
	HamiltonianMonteCarloWalk::Walk <Point, Polytope, RandomNumberGenerator, GradientFunctor, EvaluationFunctor, Solver>
	  hmc(&P, x0, grad_g, g, hmc_params);

	// Check and evaluate for all samples breaks when variance > 1, i.e. alpha > 1

	// for (int i = 1; i <= m; i++) { 					// for exact m outer loop runs
	while (alpha < 1) {									// for making the loop exit at alpha > 1

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
		std::cerr << "After i_th round | alpha = " << alpha << " | alpha_prev = " << alpha_prev << " | W_current = " << W_current << " | exp(log_W) = " << exp(log_W) << std::endl;
		alpha_prev = alpha;
	}

	return exp(log_W);    
}

template
<
	typename EvaluationFunctor,
	typename GradientFunctor,
	typename Parameters,
	typename WalkType,
	typename Polytope,
	typename Point,
	typename MT,
	typename VT,
	typename NT
>
std::pair<Point, NT> lovasz_vempala_optimize(EvaluationFunctor &g,
											GradientFunctor &grad_g,
											Parameters &params,
											Polytope &P,
											Point x0,
											NT beta = 1.0,
											NT C0 = 1,
											unsigned int walk_length = 10,
											NT delta = 1,
											NT epsilon = 0.1)
{
	unsigned int n = P.dimension();
	NT B = (NT) n * log(2 / beta);
	unsigned int m = (unsigned int)	ceil(sqrt(n) * log(2 * B * (n + log(1 / delta)) / epsilon ));
	unsigned int k = (unsigned int) ceil(C0 * n * pow(log(n), 5));

	std::cerr << "n = " << n << " m = " << m << " k = " << k << " B = " << B <<  std::endl;

	RandomNumberGenerator rng(1);
	
	// Burning samples for proper mixing
	typename WalkType::template Walk <Polytope, RandomNumberGenerator> walk(P, x0, rng);	  
	for (int i = 1; i <= k; i++) {
		walk.apply(P, x0, walk_length, rng);
	}

	// Initialize HMC walks using EvaluationFunctor and GradientFunctor
	typedef LeapfrogODESolver<Point, NT, Polytope, GradientFunctor> Solver;
	HamiltonianMonteCarloWalk::parameters <NT, GradientFunctor> hmc_params(grad_g, n);
	HamiltonianMonteCarloWalk::Walk<Point, Polytope, RandomNumberGenerator, GradientFunctor, EvaluationFunctor, Solver>
	  hmc(&P, x0, grad_g, g, hmc_params);

	
	NT alpha = (NT) 1 / B;
	Point max_point = x0;
	NT calculated_value = exp(-g(hmc.x));
	NT max_value = calculated_value;
	MT hmc_points(n, k); // n(=dimension) rows * k(=no. of points) columns 

	// Check and evaluate for all samples breaks when variance > 1, i.e. alpha > 1

	// for (int i = 1; i <= m; i++) { 					// for exact m outer loop runs
	while (alpha < 1) {									// for making the loop exit at alpha > 1

		// alpha *= (1 + 1 / sqrt(n));
		// params.set_temperature(alpha);		

		for (int j = 1; j <= k ; j++) {

			hmc.apply(rng, walk_length);
			hmc_points.col(j-1) = hmc.x.getCoefficients();

			calculated_value = exp(-g(hmc.x));
			if (calculated_value > max_value) {
				max_point = hmc.x;
				max_value = calculated_value;
			}
			
		}
	
		// MT manipulation for covariance matrix
		MT mean = hmc_points.rowwise().mean();
		for (int j = 1; j <= k ; j++) hmc_points.col(j-1) -= mean;
		
		Eigen::BDCSVD<MT> svd(hmc_points, Eigen::ComputeThinU);
		MT sigma = svd.singularValues().array().matrix().asDiagonal();

		MT covariance_matrix_inv = (1/(NT)k * svd.matrixU() * sigma * sigma * svd.matrixU().transpose()).inverse();
		VT eigen_values = covariance_matrix_inv.eigenvalues().real().array();

		NT L = eigen_values.maxCoeff();
		NT M = eigen_values.minCoeff();
		NT kappa = L/M;

		// params.update_covariance_matrix(covariance_matrix_inv, L , M, kappa)

	}

	std::pair<Point, NT> values(max_point, max_value); 
	std::cerr << "Max point check = " ; values.first.print();
	std::cerr << "Max value check = " << values.second << std::endl;

	return values;  // (Point)values.first = X s.t. max_f = max( exp(-g(X)) ), (NT)values.second = max_f

}

#endif
