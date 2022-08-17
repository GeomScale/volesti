// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer
// of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <typeinfo>
#include <unistd.h>
#include <vector>

#include "Eigen/Eigen"

#include "generators/known_polytope_generators.h"
#include "ode_solvers/ode_solvers.hpp"
#include "random.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_int.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_sequence_of_balls.hpp"

struct CustomFunctor {

  // Custom density with neg log prob equal to || x ||^2 + 1^T x
  template <typename NT> struct parameters {
    unsigned int order;
    NT L;     // Lipschitz constant for gradient
    NT m;     // Strong convexity constant
    NT kappa; // Condition number

    parameters() : order(2), L(4), m(4), kappa(1){};
  };

  template <typename Point> struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> &params;

    GradientFunctor(parameters<NT> &params_) : params(params_){};

    // The index i represents the state vector index
    Point operator()(unsigned int const &i, pts const &xs, NT const &t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * Point::all_ones(xs[0].dimension());
        y = y + (-4.0) * xs[0];
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }
  };

  template <typename Point> struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT> &params;

    FunctionFunctor(parameters<NT> &params_) : params(params_){};

    // The index i represents the state vector index
    NT operator()(Point const &x) const { return 2 * x.dot(x) + x.sum(); }
  };
};

template <typename NT> void run_main() {
  typedef Cartesian<NT> Kernel;
  typedef typename Kernel::Point Point;
  typedef std::vector<Point> pts;
  typedef HPolytope<Point> Hpolytope;
  typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
  typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
  typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
  typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor>
      Solver;

  CustomFunctor::parameters<NT> params;

  NegativeGradientFunctor F(params);
  NegativeLogprobFunctor f(params);

  RandomNumberGenerator rng(1);
  unsigned int dim = 2;

  HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(
      F, dim);

  Hpolytope P = generate_cube<Hpolytope>(dim, false);

  Point x0 = -0.25 * Point::all_ones(dim);

  HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator,
                                  NegativeGradientFunctor,
                                  NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

  int n_samples = 500; // Half will be burned

  for (int i = 0; i < n_samples; i++) {
    hmc.apply(rng, 3);
    if (i > n_samples / 2)
      std::cout << hmc.x.getCoefficients().transpose() << std::endl;
    // std::cout << hmc.solver->eta << std::endl;
  }

  std::cerr << "Step size (final): " << hmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << hmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Log-prob: "
            << exp(hmc.average_acceptance_log_prob) << std::endl;
}

int main() {
  run_main<double>();
  return 0;
}
