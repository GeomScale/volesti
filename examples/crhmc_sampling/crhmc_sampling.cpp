// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with
// Riemannian Hamiltonian
// Monte Carlo in a Constrained Space"
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <typeinfo>
#include <unistd.h>
#include <vector>

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/known_polytope_generators.h"
#include "ode_solvers/ode_solvers.hpp"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "random.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_int.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
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
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using CrhmcProblem = crhmc_problem<Point>;
  using Input = crhmc_input<MT, NT>;
  typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
  typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
  typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
  typedef ImplicitMidpointODESolver<Point, NT, CrhmcProblem,
                                    NegativeGradientFunctor>
      Solver;
  using Opts = opts<NT>;

  CustomFunctor::parameters<NT> params;

  NegativeGradientFunctor F(params);
  NegativeLogprobFunctor f(params);

  RandomNumberGenerator rng(1);
  unsigned int dim = 2;
  Opts options;
  CRHMCWalk::parameters<NT, NegativeGradientFunctor> crhmc_params(F, dim,
                                                                  options);

  Input input = Input(dim);
  input.lb = -VT::Ones(dim);
  input.ub = VT::Ones(dim);
  CrhmcProblem P = CrhmcProblem(input);
  Point x0 = Point(P.center);

  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator,
                  NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      crhmc(P, x0, F, f, crhmc_params);

  int n_samples = 50000; // Half will be burned
  for (int i = 0; i < n_samples; i++) {
    crhmc.apply(rng, 3);
    if (i > n_samples / 2)
      std::cout << crhmc.x.getCoefficients().transpose() << std::endl;
    // std::cout << hmc.solver->eta << std::endl;
  }

  std::cerr << "Step size (final): " << crhmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << crhmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Log-prob: "
            << exp(crhmc.average_acceptance_log_prob) << std::endl;
}

int main() {
  run_main<double>();
  return 0;
}
