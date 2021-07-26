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
typedef std::vector<Point> Points;
typedef HPolytope<Point> HPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
typedef LinearProgramFunctor::GradientFunctor<Point> NegativeGradientFunctor;
typedef LinearProgramFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
typedef OptimizationFunctor::GradientFunctor
  <Point, NegativeLogprobFunctor, NegativeGradientFunctor> NegativeGradientOptimizationFunctor;
typedef OptimizationFunctor::FunctionFunctor
  <Point, NegativeLogprobFunctor, NegativeGradientFunctor> NegativeLogprobOptimizationFunctor;

typedef const unsigned int Uint;  // Positive constant value for no of samples & dimensions
enum volumetype { CB ,CG ,SOB }; // Volume type for polytope

template
<
    typename EvaluationFunctor,
    typename GradientFunctor,
    typename WalkType,
    typename Polytope = HPOLYTOPE,
    typename RNG = RandomNumberGenerator,
    typename MT,
    typename VT,
    typename NT
>
NT lovasz_vempala_integrate(EvaluationFunctor &g,
                            GradientFunctor &grad_g,
                            Polytope &P,
                            Point x0,
                            NT B,
                            unsigned int walk_length = 10,
                            NT epsilon = 0.1)
{
    unsigned int n = P.dimension();
    unsigned int m = (unsigned int) ceil(sqrt(n) * log(B));
    unsigned int k = (unsigned int) ceil(512 / pow(epsilon,2) * sqrt(n) * log(B));
    
    NT volume = volume_sequence_of_balls <BallWalk, RNG, HPOLYTOPE> (P, epsilon, walk_length);
    NT a = (NT) 1 / B;
    NT alpha = 0, alpha_prev = 0;
    NT log_W = log(volume);
    NT W_current = (NT) 0;

    // MT points = uniform sample points in ` n rows * k columns ` ( k is the number of points )
    RNG rng(1);
    MT points(n,k); // n rows(dimension) * k columns(number of points)
    typename WalkType::template Walk <Polytope, RNG> walk(P, x0, rng);
     
    for (int i = 0; i < k; i++) {
        walk.apply(P, x0, walk_length, rng);
        points.col(i) = x0.getCoefficients();
    }

    // Initialize ReHMC using OptimizationFunctor
    LinearProgramFunctor::parameters<NT, Point> lp_params(x0);
    NegativeGradientFunctor F_lp(lp_params);
    NegativeLogprobFunctor f_lp(lp_params);

    // Declare optimization oracles using g and grad_g
    typedef LeapfrogODESolver<Point, NT, Polytope,  NegativeGradientOptimizationFunctor> Solver;

    OptimizationFunctor::parameters
      <NT, NegativeLogprobFunctor, NegativeGradientFunctor>
      opt_params(1, x0.dimension(), f_lp, F_lp);

    NegativeLogprobOptimizationFunctor f(opt_params);
    NegativeGradientOptimizationFunctor F(opt_params);

    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientOptimizationFunctor> hmc_params(F, n);
    
    HamiltonianMonteCarloWalk::Walk
      <Point, Polytope, RandomNumberGenerator, NegativeGradientOptimizationFunctor, NegativeLogprobOptimizationFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    // Check and evaluate for all samples breaks when variance > 1, i.e. a > 1
    int i = 1;
    while ( i++ <= m || a <= 1 ) {

        a *= (1 + 1 / sqrt(n)); // variance sequence algorithm stops when variance > 1
        alpha = a;
        // g_var.set_variance(a);
        // grad_g_var.set_variance(a);
        W_current = 0;

        for (unsigned int j = 1; j <= k ; j++) {

            hmc.apply(rng, walk_length);
            W_current += exp(-g(hmc.x) * (alpha - alpha_prev));

        }

        W_current /= k;
        log_W += log(W_current);
        alpha_prev = a;
    }

    return exp(log_W);
    
}

// include/ode_solvers/oracle_functors 
// To see how to define oracle_functors and gradient functors
// https://github1s.com/GeomScale/volume_approximation/blob/develop/include/ode_solvers/oracle_functors.hpp

// test/logconcave_sampling_test.cpp
// https://github1s.com/GeomScale/volume_approximation/blob/develop/test/logconcave_sampling_test.cpp

// HMC examples/logconcave
// https://github.com/GeomScale/volume_approximation/tree/develop/examples/logconcave

/*
template 
<
    typename EvaluationFunctor,
    typename NT
>
struct VarianceEvaluationFunctor {
    EvaluationFunctor &g;
    NT a;

    VarianceEvaluationFunctor(EvaluationFunctor g_, NT a_ ) : g(g_) , a(a_) {} ;

	NT operator()(Point x) {// see the oracle_functors on how to define this operator 
		return a * g(x);
	}

    void set_variance (NT a_) {
        a = a_;
    }
};

template
<
    typename GradientFunctor,
    typename Point
>
struct VarianceGradientFunctor {

    GradientFunctor &g;
    NT a;

    VarianceGradientFunctor(GradientFunctor g_, NT a_ ) : g(g_), a(a_) {} ; 

    Point operator()(Point x) {// see the oracle_functors on how to define this operator 
		return a * g(x);
	}

    void set_variance(NT a_) {
        a = a_;
    }

};
*/

/*
    // Declare LinearProgramFunctor as oracles for OptimizationFunctor
    // LinearProgramFunctor::parameters<NT, Point> lp_params(x0);
    // NegativeGradientFunctor F_lp(lp_params);
    // NegativeLogprobFunctor f_lp(lp_params);

    OptimizationFunctor::parameters <NT, VarianceGradientFunctor> opt_params(1, x0.dimension(), g_var, g_grad_var);
    
    VarianceEvaluationFunctor <EvaluationFunctor, NT> g_var(opt_params); //f
    VarianceGradientFunctor <GradientFunctor, Point> g_grad_var(opt_params); //F    

    HamiltonianMonteCarloWalk::parameters
      <NT, VarianceGradientFunctor> hmc_params(g_grad_var, n);

    HamiltonianMonteCarloWalk::Walk <Point, Polytope, RNG, VarianceGradientFunctor, Solver>
      hmc(&P, x0, g_grad_var, g_var, hmc_params);
*/