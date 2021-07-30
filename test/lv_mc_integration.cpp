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
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "misc.h"

typedef double NT;

template <typename NT>
void test_values (NT computed){ //, NT expected, NT exact) {
	std::cout << "Computed integration value = " << computed << std::endl;
	// std::cout << "Expected integration value = " << expected << std::endl;
	// std::cout << "Exact integration value = " << exact << std::endl;
	// std::cout << "Relative error (expected) = " << std::abs((computed - expected)/expected) << std::endl;
	// std::cout << "Relative error (exact) = " << std::abs((computed - exact)/exact) << std::endl ;
	// CHECK(((std::abs((computed - expected)/expected) < 0.00001) || (std::abs((computed - exact)/exact) < 0.2)));
}

struct CustomFunctor {

  // Custom density with neg log prob equal to || x ||^2 + 1^T x
  template 
  <
      typename NT
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number

    parameters() : order(2), L(2), m(2), kappa(1) {};

    parameters(unsigned int order_) :
      order(order),
      L(2),
      m(2),
      kappa(1)
    {}
  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> params;

    GradientFunctor() {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * Point::all_ones(xs[0].dimension());
        y = y + (-2.0) * xs[0];
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT> params;

    FunctionFunctor() {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      return x.dot(x) + x.sum();
    }

  };

};

template <typename NT>
void call_cubes_test_lovasz_vempala_integrate() { // or inside the previous test function

	typedef Cartesian<NT> Kernel;
	typedef typename Kernel::Point Point;
	typedef HPolytope<Point> HPOLYTOPE;
	typedef boost::mt19937 RNGType;
	typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;

	typedef CustomFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef CustomFunctor::GradientFunctor <Point> GradientFunctor;
	typedef OptimizationFunctor::GradientFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeGradientOptimizationFunctor;
	typedef OptimizationFunctor::FunctionFunctor
	<Point, EvaluationFunctor, GradientFunctor> NegativeLogprobOptimizationFunctor;
	typedef OptimizationFunctor::parameters<NT, EvaluationFunctor, GradientFunctor> OptimizationParameters;

	GradientFunctor grad_g;
	EvaluationFunctor g;
	HPOLYTOPE HP = generate_cube <HPOLYTOPE> (2, false);
	// HP.print();

	OptimizationParameters opt_params(1, HP.dimension(), g, grad_g);

	NegativeLogprobOptimizationFunctor f(opt_params);
	NegativeGradientOptimizationFunctor grad_f(opt_params);

	std::vector<NT> Origin{0,0}; Point x1(2, Origin);
	std::vector<NT> Corner{1,1}; Point x2(2, Corner);

	NT B = log( exp(-g(x1)) / exp(-g(x2)) ) ; // 2*n + 2*log(1/0.1) + n*log( 1 / beta);;
	NT integral_value;
	unsigned int n = HP.dimension();;
	std::pair <Point, NT> inner_ball = HP.ComputeInnerBall();;
	Point x0 = inner_ball.first;

	integral_value = lovasz_vempala_integrate <NegativeLogprobOptimizationFunctor, NegativeGradientOptimizationFunctor, OptimizationParameters, HPOLYTOPE, Point, NT>
	  (f, grad_f, opt_params, HP, x0, B, 10, 0.1);
	
	test_values(integral_value);

}

TEST_CASE("cubes") {
    call_cubes_test_lovasz_vempala_integrate<double>();
}
