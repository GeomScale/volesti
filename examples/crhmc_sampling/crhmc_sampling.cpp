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
#include "diagnostics/multivariate_psrf.hpp"
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
  // Custom Gaussian density
  template <typename NT> struct parameters {
    unsigned int order;
    NT L;     // Lipschitz constant for gradient
    NT m;     // Strong convexity constant
    NT kappa; // Condition number
    parameters() : L(4), m(4), kappa(1){};
  };

  template <typename Point> struct Grad {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;
    parameters<NT> &params;
    Grad(parameters<NT> &params_) : params(params_){};
    Point operator()(Point const &x) const {
      Point y = -(1.0) * x;
      return y;
    }
  };
  template <typename Point> struct Hess {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> &params;
    Hess(parameters<NT> &params_) : params(params_){};
    Point operator()(Point const &x) const {
      return Point::all_ones(x.dimension());
    }
  };

  template <typename Point> struct Func {
    typedef typename Point::FT NT;
    parameters<NT> &params;
    Func(parameters<NT> &params_) : params(params_){};
    NT operator()(Point const &x) const { return 0.5 * x.dot(x); }
  };
};

template <typename NT> void run_main(int n_samples = 500, int n_burns = -1) {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using pts = std::vector<Point>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;
  using Func = CustomFunctor::Func<Point>;
  using Grad = CustomFunctor::Grad<Point>;
  using Hess = CustomFunctor::Hess<Point>;
  using Input = crhmc_input<MT, Point, Func, Grad>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Solver = ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad>;
  using Opts = opts<NT>;

  CustomFunctor::parameters<NT> params;

  Func f(params);
  Grad g(params);
  Hess h(params);
  if (n_burns == -1) {
    n_burns = n_samples / 2;
  }
  RandomNumberGenerator rng(1);
  unsigned int dim = 2;
  Opts options;
  CRHMCWalk::parameters<NT, Grad> crhmc_params(g, dim, options);
  MT A = MT::Ones(5, dim);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = 10 * VT::Ones(5, 1);
  Input input = Input(dim, f, g);
  input.Aineq = A;
  input.bineq = b;
  CrhmcProblem P = CrhmcProblem(input);
  P.print();
  Point x0 = Point(P.center);
  crhmc_params.eta = 0.2;
  crhmc_params.momentum = 0.8;
  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator, Grad, Func,
                  Solver>
      crhmc(P, x0, g, f, crhmc_params);
  MT samples = MT(dim, n_samples - n_burns);
  int j = 0;
  for (int i = 0; i < n_samples; i++) {
    crhmc.apply(rng, 30, true);
    if (i > n_burns) {
      VT sample = P.T * crhmc.x.getCoefficients() + P.y;
      samples.col(j) = VT(sample);
      j++;
      std::cout << sample.transpose() << std::endl;
    }
  }

  std::cerr << "Step size (final): " << crhmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << crhmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Probability: "
            << exp(crhmc.average_acceptance_log_prob) << std::endl;
  std::cerr << "PSRF: " << multivariate_psrf<NT, VT, MT>(samples) << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc == 1)
    run_main<double>();
  else if (argc == 2)
    run_main<double>(atoi(argv[1]));
  else if (argc == 3)
    run_main<double>(atoi(argv[1]), atoi(argv[2]));

  return 0;
}
