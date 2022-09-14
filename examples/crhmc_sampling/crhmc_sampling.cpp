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
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "diagnostics/multivariate_psrf.hpp"
#include "generators/known_polytope_generators.h"
#include "misc/misc.h"
#include "ode_solvers/ode_solvers.hpp"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "random.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_int.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <typeinfo>
#include <unistd.h>
#include <vector>

template <typename NT, int simdLen = 1>
void run_main(int n_samples = 10000, int n_burns = -1, int dimension = 2,
              int walk_length = 1, int burn_steps = 1) {
  std::cerr << "Using simdLen= " << simdLen << "\n";
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using pts = std::vector<Point>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;
  using Func = GaussianFunctor::FunctionFunctor<Point>;
  using Grad = GaussianFunctor::GradientFunctor<Point>;
  using Hess = GaussianFunctor::HessianFunctor<Point>;
  using func_params = GaussianFunctor::parameters<NT, Point>;
  using Input = crhmc_input<MT, Point, Func, Grad, Hess>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Solver =
      ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad, simdLen>;
  using Opts = opts<NT>;
  RandomNumberGenerator rng(1);
  if (n_burns == -1) {
    n_burns = n_samples / 2;
  }
  func_params params = func_params(Point(dimension), 4, 1);
  Func f(params);
  Grad g(params);
  Hess h(params);
  Opts options;
  options.simdLen = simdLen;
  Input input = Input(dimension, f, g, h);
  input.lb = -VT::Ones(dimension);
  input.ub = VT::Ones(dimension);
  CrhmcProblem P = CrhmcProblem(input, options);
  P.print();
  Point x0 = Point(P.center);
  CRHMCWalk::parameters<NT, Grad> crhmc_params(g, dimension, options);
  crhmc_params.eta = 0.2;
  crhmc_params.momentum = 0.8;
  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator, Grad, Func,
                  Solver>
      crhmc(P, x0, g, f, crhmc_params);
  int max_actual_draws = n_samples - n_burns;
  MT samples = MT(dimension, max_actual_draws);
#ifdef TIME_KEEPING
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::system_clock::now();
#endif
  for (int i = 0; i < n_samples; i++) {
    if (i % 1000 == 0) {
      std::cerr << i << " out of " << n_samples << "\n";
    }
    for (int k = 0; k < burn_steps; k++) {
      crhmc.apply(rng, walk_length, true);
    }
    if (i >= n_burns) {
      MT sample = crhmc.getPoints();
      int j = i - n_burns;
      if (j * simdLen + simdLen - 1 < max_actual_draws) {
        samples(Eigen::all,
                Eigen::seq(j * simdLen, j * simdLen + simdLen - 1)) = sample;
      } else {
        samples(Eigen::all, Eigen::seq(j * simdLen, max_actual_draws - 1)) =
            sample(Eigen::all,
                   Eigen::seq(0, max_actual_draws - 1 - simdLen * j));
      }
    }
  }
  std::cerr << "\n";
#ifdef TIME_KEEPING
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> total_time = end - start;
  std::cerr << "Total time: " << total_time.count() << "\n";
  crhmc.print_timing_information(std::cerr);
#endif

  std::cerr << "Step size (final): " << crhmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << crhmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Probability: "
            << crhmc.average_acceptance_prob << std::endl;
  std::cerr << "PSRF: " << multivariate_psrf<NT, VT, MT>(samples) << std::endl;
  std::cerr << "Writing samples in a file \n";
#ifdef TIME_KEEPING
  start = std::chrono::system_clock::now();
#endif
  std::cout << samples.transpose() << std::endl;
#ifdef TIME_KEEPING
  end = std::chrono::system_clock::now();
  total_time = end - start;
  std::cerr << "Time for writing the file: " << total_time.count() << "\n";
#endif
}

int main(int argc, char *argv[]) {
  std::cerr
      << "Example Usage: ./crhmc_sampling [n_samples] [initial_burns] "
         "[dimension] [ode_steps] [steps_bettween_samples] [simdLen 1 or 4]\n";
  std::cerr << "Example Usage: ./crhmc_sampling 10000 5000 "
               "2 1 1 4 > samples.txt\n";
  if (argc == 1) {
    run_main<double>();
  } else if (argc == 2) {
    run_main<double>(atoi(argv[1]));
  } else if (argc == 3) {
    run_main<double>(atoi(argv[1]), atoi(argv[2]));
  } else if (argc == 4) {
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
  } else if (argc == 5) {
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
                     atoi(argv[4]));
  } else if (argc == 6) {
    run_main<double>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
                     atoi(argv[5]));
  } else if (argc == 7) {
    int simdLen = atoi(argv[6]);
    if (simdLen == 1) {
      run_main<double, 1>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
                          atoi(argv[4]), atoi(argv[5]));
    } else if (simdLen == 4) {
      run_main<double, 4>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
                          atoi(argv[4]), atoi(argv[5]));
    }
  }
  return 0;
}
