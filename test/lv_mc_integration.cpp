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
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef std::vector<Point> Points;
typedef HPolytope<Point> HPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;

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
  template <
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

template <typename NT = NT>
void call_cubes_test_lovasz_vempala_integrate() { // or inside the previous test function

  typedef CustomFunctor::FunctionFunctor <Point> EvaluationFunctor;
	typedef CustomFunctor::GradientFunctor <Point> GradientFunctor;

  CustomFunctor::parameters<NT> params;

  GradientFunctor grad_g;
  EvaluationFunctor g;
	HPOLYTOPE HP;

	HP = generate_cube <HPOLYTOPE> (2, false);

  std::vector<NT> Origin{0, 0};
	Point x0(2, Origin);
  std::vector<NT> Corner{1, 1};
	Point x1(2, Corner);

	NT beta = 1; //log(exp(-g(x0))) / exp(-g(x1));
	unsigned int n = HP.dimension();
	NT B = log( exp(-g(x0)) / exp( -g(x1)) ) ; // 2*n + 2*log(1/0.1) + n*log( 1 / beta);

	NT integral_value = lovasz_vempala_integrate <EvaluationFunctor, GradientFunctor,BilliardWalk, HPOLYTOPE, Point, NT>
	  (g, grad_g, HP, x0, B, 10, 0.1);
	
  test_values(integral_value);

}

TEST_CASE("cubes"){
    call_cubes_test_lovasz_vempala_integrate();
}

/*
    e.g. use the function f(x) = exp( -g(x))
    g(x) returns ||x||^2
    f(x) return exp(-g(x))

    NT result1 = simple_mc_integrate <> ( ... )
    NT result2 = lovasz_vempala_integrate <> ( ... )

    Polytope K ; // eg use a cube
    Max_f = exp( -g(0) );
    Min_f = exp( -g(point at the corner) );
    B = ...
    Point x0 = K.inner_ball();

    CHECK ( relation_error(result1, result2) is small )
*/
