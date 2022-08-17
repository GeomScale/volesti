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

    parameters() : order(1), L(4), m(4), kappa(1){};
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
        // y=Point(xs[0].dimension())-xs[0];
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }
    Point operator()(Point const &x) const {
      Point y = (-1.0) * Point::all_ones(x.dimension());
      y = y + (-4.0) * x;
      return y;
    }
  };
  template <typename Point> struct HessianFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> &params;
    HessianFunctor(parameters<NT> &params_) : params(params_){};
    Point operator()(Point const &x) const { return 2 * Point(x.dimension()); }
  };

  template <typename Point> struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT> &params;

    FunctionFunctor(parameters<NT> &params_) : params(params_){};
    // The index i represents the state vector index
    NT operator()(Point const &x) const { return 2 * x.dot(x) + x.sum(); }
    //  NT operator()(Point const &x) const { return (0.5) * x.dot(x); }
  };
};

template <typename NT> void run_main(int n_samples = 500) {
  typedef Cartesian<NT> Kernel;
  typedef typename Kernel::Point Point;
  typedef std::vector<Point> pts;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
  typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
  typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
  typedef CustomFunctor::HessianFunctor<Point> HessianFunctor;
  using Input =
      crhmc_input<MT, Point, NegativeLogprobFunctor, NegativeGradientFunctor>;

  using CrhmcProblem = crhmc_problem<Point, Input>;
  typedef ImplicitMidpointODESolver<Point, NT, CrhmcProblem,
                                    NegativeGradientFunctor>
      Solver;
  using Opts = opts<NT>;

  CustomFunctor::parameters<NT> params;

  NegativeGradientFunctor F(params);
  NegativeLogprobFunctor f(params);
  HessianFunctor h(params);
  RandomNumberGenerator rng(1);
  unsigned int dim = 2;
  Opts options;
  CRHMCWalk::parameters<NT, NegativeGradientFunctor> crhmc_params(F, dim,
                                                                  options);

  MT A = MT::Ones(5, dim);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = 10 * VT::Ones(5, 1);
  Input input = Input(dim, f, F);
  input.Aineq = A;
  input.bineq = b;
  CrhmcProblem P = CrhmcProblem(input);
  P.print();
  Point x0 = Point(P.center);
  crhmc_params.eta = 0.2;
  crhmc_params.momentum = 0.8;
  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator,
                  NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      crhmc(P, x0, F, f, crhmc_params);

  for (int i = 0; i < n_samples; i++) {
    crhmc.apply(rng, 30, true);
    if (i > n_samples / 2)
      std::cout << (P.T * crhmc.x.getCoefficients() + P.y).transpose()
                << std::endl;
  }

  std::cerr << "Step size (final): " << crhmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << crhmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Log-prob: "
            << exp(crhmc.average_acceptance_log_prob) << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc == 1)
    run_main<double>();
  else
    run_main<double>(atoi(argv[1]));

  return 0;
}
