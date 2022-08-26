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
#include "misc/misc.h"
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
    NT var = 18;
    parameters() : L(4), m(4), kappa(1){};
  };

  template <typename Point> struct Grad {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;
    parameters<NT> &params;
    Grad(parameters<NT> &params_) : params(params_){};
    Point operator()(Point const &x) const {
      Point y = -(1.0 / params.var) * x;
      return y;
    }
  };
  template <typename Point> struct Hess {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> &params;
    Hess(parameters<NT> &params_) : params(params_){};
    Point operator()(Point const &x) const {
      return (1.0 / params.var) * Point::all_ones(x.dimension());
    }
  };

  template <typename Point> struct Func {
    typedef typename Point::FT NT;
    parameters<NT> &params;
    Func(parameters<NT> &params_) : params(params_){};
    NT operator()(Point const &x) const {
      return (1.0 / params.var) * 0.5 * x.dot(x);
    }
  };
};

template <typename NT>
void run_main(int n_samples = 10000, int n_burns = -1, int dimension = 2,
              int walk_length = 1, int burn_steps = 1, int dynamic_weight = 0,
              int dynamic_regularization = 0, int dynamic_step_size = 0) {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using pts = std::vector<Point>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;
  using Func = CustomFunctor::Func<Point>;
  using Grad = CustomFunctor::Grad<Point>;
  using Hess = CustomFunctor::Hess<Point>;
  // using Input = crhmc_input<MT, Point>;
  using Input = crhmc_input<MT, Point, Func, Grad, Hess>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Solver = ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad>;
  using Opts = opts<NT>;
  using Hpolytope = HPolytope<Point>;

  CustomFunctor::parameters<NT> params;

  Func f(params);
  Grad g(params);
  Hess h(params);
  if (n_burns == -1) {
    n_burns = n_samples / 2;
  }
  RandomNumberGenerator rng(1);
  unsigned int dim = dimension;
  /*
    std::cerr << "Reading input from file..." << std::endl;
    std::ifstream inp;
    std::vector<std::vector<NT>> Pin;
    inp.open("../../test/metabolic_full_dim/polytope_e_coli.ine",
             std::ifstream::in);
    read_pointset(inp, Pin);
    Hpolytope Polytope(Pin);
    dim = Polytope.dimension();
  */
  Opts options;
  if (dynamic_weight == 1) {

    options.DynamicWeight = true;
  } else {
    options.DynamicWeight = false;
  }
  if (dynamic_regularization == 1) {

    options.DynamicRegularizer = true;
  } else {
    options.DynamicRegularizer = false;
  }
  if (dynamic_step_size == 1) {
    options.DynamicStepSize = true;
  } else {
    options.DynamicStepSize = false;
  }
  CRHMCWalk::parameters<NT, Grad> crhmc_params(g, dim, options);
  MT A = MT::Ones(5, dim);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = 10 * VT::Ones(5, 1);
  // Hpolytope Polytope = generate_simplex<Hpolytope>(dim, false);
  // MT A = Polytope.get_mat();
  // VT b = Polytope.get_vec();
  // std::cerr<<"A.rows============== " << A.rows()<<"\n";
  // std::cerr<<"A=\n"<<A<<"\n";
  // std::cerr<<"b=\n"<<b.transpose()<<"\n";

  Input input = Input(dim, f, g, h);
  input.Aineq = A;
  input.bineq = b;
  // input.lb = -VT::Ones(dim);
  // input.ub = VT::Ones(dim);
  CrhmcProblem P = CrhmcProblem(input, options);
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
    if (i % 1000 == 0) {
      std::cerr << i << " out of " << n_samples << "\n";
    }
    for (int k = 0; k < burn_steps; k++) {
      crhmc.apply(rng, walk_length, true);
    }
    if (i >= n_burns) {
      VT sample = crhmc.getPoint().getCoefficients();
      samples.col(j) = VT(sample);
      j++;
      std::cout << sample.transpose() << std::endl;
    }
  }
  std::cerr << "\n";

  std::cerr << "Step size (final): " << crhmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << crhmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Probability: "
            << exp(crhmc.average_acceptance_log_prob) << std::endl;
  std::cerr << "PSRF: " << multivariate_psrf<NT, VT, MT>(samples) << std::endl;
}

int main(int argc, char *argv[]) {
  std::cerr << "Example Usage: ./crhmc_sampling n_sample initial_burns "
               "dimension ode_steps steps_bettween_samples\n";
  std::cerr << "Example Usage: ./crhmc_sampling 10000 5000 "
               "2 1 1\n";
  if (argc == 1)
    run_main<double>();
  else if (argc == 2)
    run_main<double>(atoi(argv[1]));
  else if (argc == 3)
    run_main<double>(atoi(argv[1]), atoi(argv[2]));
  else if (argc == 4)
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
  else if (argc == 5)
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
                     atoi(argv[4]));
  else if (argc == 6)
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                     atoi(argv[5]));
  else if (argc == 7)
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                     atoi(argv[5]), atoi(argv[6]));
  else if (argc == 8)
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                     atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
  else if (argc == 9)
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                     atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
                     atoi(argv[8]));
  return 0;
}
