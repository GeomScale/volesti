// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Optimization algorithm used here : https://www.cc.gatech.edu/~vempala/acg/www/www/papers/integration.pdf (Page-8)

#ifndef LOVASZ_VEMPALA_OPTIMIZATION_HPP
#define LOVASZ_VEMPALA_OPTIMIZATION_HPP

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
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> HPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;

template
<
	typename MT,
	typename NT
>
MT covariance_matrix_calculator(MT &matrix) {

	unsigned int k = matrix.cols();
	MT mean = matrix.rowwise().mean();
	for (int j = 1; j <= k; j++) matrix.col(j-1) -= mean;
	
	Eigen::BDCSVD<MT> svd(matrix, Eigen::ComputeThinU);
	MT sigma = svd.singularValues().array().matrix().asDiagonal();

	MT covariance_matrix = (NT)1/(NT)k * svd.matrixU() * sigma * sigma * svd.matrixU().transpose();
	return covariance_matrix;

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
	MT points_matrix(n, k); // n(=dimension) rows * k(=no. of points) columns 
	
	// Burning samples for proper mixing
	typename WalkType::template Walk <Polytope, RandomNumberGenerator> walk(P, x0, rng);	  
	for (int i = 1; i <= k; i++) {
		walk.apply(P, x0, walk_length, rng);
		points_matrix.col(i-1) = x0.getCoefficients();
	}

	MT covariance_matrix = covariance_matrix_calculator<MT, NT>(points_matrix);
	params.update_covariance_matrix(covariance_matrix);

	// Initialize HMC walks using EvaluationFunctor and GradientFunctor
	typedef LeapfrogODESolver<Point, NT, Polytope, GradientFunctor> Solver;
	HamiltonianMonteCarloWalk::parameters <NT, GradientFunctor> hmc_params(grad_g, n);
	HamiltonianMonteCarloWalk::Walk<Point, Polytope, RandomNumberGenerator, GradientFunctor, EvaluationFunctor, Solver>
	  hmc(&P, x0, grad_g, g, hmc_params);

	
	NT alpha = (NT) 1 / B;
	Point max_point = x0;
	NT calculated_value = exp(-g(hmc.x));
	NT max_value = calculated_value;

	// Check and evaluate for all samples breaks when variance > 1, i.e. alpha > 1

	for (int i = 1; i <= m; i++) { 					// for exact m outer loop runs
	// while (alpha < 1) {							// for making the loop exit at alpha > 1

		// alpha *= (1 + 1 / sqrt(n));
		// params.set_temperature(alpha);		

		for (int j = 1; j <= k ; j++) {

			hmc.apply(rng, walk_length);
			points_matrix.col(j-1) = hmc.x.getCoefficients();

			calculated_value = exp(-g(hmc.x));
			if (calculated_value > max_value) {
				max_point = hmc.x;
				max_value = calculated_value;
			}
			
		}

		covariance_matrix = covariance_matrix_calculator<MT, NT>(points_matrix);
		params.update_covariance_matrix(covariance_matrix);
	}

	std::pair<Point, NT> values(max_point, max_value); 
	std::cerr << "Max point check = " ; values.first.print();
	std::cerr << "Max value check = " << values.second << std::endl;
	
	/*	
		(Point)values.first = X s.t. max_f = max( exp(-g(X)) )
		(NT)values.second = max_f
	*/
	return values;  
}

#endif
