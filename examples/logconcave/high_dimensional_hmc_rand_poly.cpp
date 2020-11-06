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

#include "Eigen/Eigen"

#include "misc/misc.h"
#include "misc/linear_extensions.h"
#include "lp_oracles/solve_lp.h"
#include "ode_solvers/ode_solvers.hpp"
#include "diagnostics/geweke.hpp"
#include "diagnostics/multivariate_psrf.hpp"
#include "diagnostics/raftery.hpp"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"

struct CustomFunctor {

  // Custom density with neg log prob equal to c^T x
  template <
      typename NT,
      typename Point
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number
    Point x0;

    parameters(Point x0_) : order(2), L(1), m(1), kappa(1), x0(x0_) {};

  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT, Point> &params;

    GradientFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * (xs[0] - params.x0);
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

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      Point y = x - params.x0;
      return 0.5 * y.dot(y);
    }

  };

};

template <typename NT>
void run_main() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef boost::mt19937 RNGType;
    typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
    typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;
    typedef typename HPolytope<Point>::MT MT;
    typedef typename HPolytope<Point>::VT VT;


    unsigned int dim = 1000;
    Hpolytope P = random_hpoly<Hpolytope, RNGType>(dim, 3 * dim);

    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();

    // Chebyshev center
    Point x0 = inner_ball.first;

    CustomFunctor::parameters<NT, Point> params(x0);

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);



    RandomNumberGenerator rng(1);

    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);


    // In the first argument put in the address of an H-Polytope
    // for truncated sampling and NULL for untruncated
    HamiltonianMonteCarloWalk::Walk
      <Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    int n_samples = 80000;
    int n_burns = 0;

    MT samples;
    samples.resize(dim, n_samples - n_burns);

    hmc.solver->eta0 = 0.5;

    for (int i = 0; i < n_samples; i++) {
      if (i % 1000 == 0) std::cerr << ".";
      hmc.apply(rng, 3);
      if (i >= n_burns) {
          samples.col(i - n_burns) = hmc.x.getCoefficients();
          std::cout << hmc.x.getCoefficients().transpose() << std::endl;
      }
    }
    std::cerr << std::endl;


    std::cerr << "Step size (final): " << hmc.solver->eta << std::endl;
    std::cerr << "Discard Ratio: " << hmc.discard_ratio << std::endl;
    std::cerr << "Average Acceptance Probability: " << exp(hmc.average_acceptance_log_prob) << std::endl;
    std::cerr << "PSRF: " <<  multivariate_psrf<NT, VT, MT>(samples) << std::endl;
}

int main() {
  run_main<double>();
  return 0;
}
