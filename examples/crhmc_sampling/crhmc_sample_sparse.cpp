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
#include "diagnostics/diagnostics.hpp"
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
#include <unsupported/Eigen/SparseExtra>
#include <vector>
template <typename MT, typename Polytope, typename RandomNumberGenerator,
          typename WalkTypePolicy, typename NT, typename Point, typename Input,
          typename Solver, typename Opts, typename StreamType>
void sample(MT &samples, Polytope &P, RandomNumberGenerator &rng,
            const unsigned int n_samples, const unsigned int n_burns,
            Input &input, Opts &options, StreamType &stream) {

  using NegativeGradientFunctor = typename Input::Grad;
  using NegativeLogprobFunctor = typename Input::Func;
  typedef typename WalkTypePolicy::template Walk<
      Point, Polytope, RandomNumberGenerator, NegativeGradientFunctor,
      NegativeLogprobFunctor, Solver>
      walk;
  typedef
      typename WalkTypePolicy::template parameters<NT, NegativeGradientFunctor>
          walk_params;
  int dimension = input.dimension;
  int simdLen = options.simdLen;
  Point p = Point(P.center);
  walk_params params(input.df, p.dimension(), options);

  if (input.df.params.eta > 0) {
    params.eta = input.df.params.eta;
  }
  walk crhmc_walk = walk(P, p, input.df, input.f, params);
  std::cerr << "Burn-in "<<n_burns<<" draws" << std::endl;
  for (int i = 0; i < n_burns; i++) {
    if (i % 1000 == 0) {
      std::cerr << i << " out of " << n_burns << "\n";
    }
    crhmc_walk.apply(rng, 1);
  }
  int max_actual_draws = n_samples - n_burns;
  samples=MT(dimension, max_actual_draws);
  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
  std::cerr << "Sampling "<<std::ceil(max_actual_draws / simdLen)<<" draws" << std::endl;
  crhmc_walk.initialize_timers();
  start = std::chrono::system_clock::now();
  for (unsigned int i = 0; i < std::ceil(max_actual_draws / simdLen); i++) {
    if (i % 1000 == 0) {
      std::cerr << i << " out of " << std::ceil(max_actual_draws / simdLen)
                << "\n";
    }
    crhmc_walk.apply(rng, 1);
    MT sample = crhmc_walk.getPoints();
    if (i * simdLen + simdLen - 1 < max_actual_draws) {
      samples(Eigen::all, Eigen::seq(i * simdLen, i * simdLen + simdLen - 1)) =
          sample;
    } else {
      samples(Eigen::all, Eigen::seq(i * simdLen, max_actual_draws - 1)) =
          sample(Eigen::all, Eigen::seq(0, max_actual_draws - 1 - simdLen * i));
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> total_time = stop - start;
  crhmc_walk.print_timing_information(stream);
  stream << "---Total Sampling time: " << total_time.count() << "\n";
  std::cerr << "---Total Sampling time: " << total_time.count() << "\n";
  stream << "Number of non Zeros: " << P.nnz() << std::endl;
  stream << "Step size (final): " << crhmc_walk.solver->eta << std::endl;
  stream << "Discard Ratio: " << crhmc_walk.discard_ratio << std::endl;
  stream << "Average Acceptance Probability: "
         << crhmc_walk.average_acceptance_prob << std::endl;
  delete crhmc_walk.module_update;
}
inline bool exists_check(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}
template <typename NT, typename VT, typename MT>
NT max_interval_psrf(MT &samples) {
  NT max_psrf = NT(0);
  VT intv_psrf = interval_psrf<VT, NT, MT>(samples);
  unsigned int d = intv_psrf.rows();
  for (unsigned int i = 0; i < d; i++) {
    if (intv_psrf(i) > max_psrf)
      max_psrf = intv_psrf(i);
  }
  return max_psrf;
}
template <typename MT,typename VT,typename NT, typename StreamType>
void diagnose(MT &samples, StreamType &stream) {
  unsigned int min_ess = 0;
  print_diagnostics<NT, VT, MT>(samples, min_ess, stream);
  max_psrf = max_interval_psrf<NT, VT, MT>(samples);
  stream << "max_psrf: " << max_psrf << std::endl;
  stream << "min ess " << min_ess << std::endl;
}
using NT = double;
using Kernel = Cartesian<NT>;
using Point = typename Kernel::Point;
using Func = GaussianFunctor::FunctionFunctor<Point>;
using Grad = GaussianFunctor::GradientFunctor<Point>;
using Hess = GaussianFunctor::HessianFunctor<Point>;
using func_params = GaussianFunctor::parameters<NT, Point>;
using Input = crhmc_input<SpMat, Point, Func, Grad, Hess>;
using CrhmcProblem = crhmc_problem<Point, Input>;
using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
using triplet = Eigen::Triplet<NT>;
using Opts = opts<NT>;
using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;
/*Problem on the form X=[A|b] bounds=[lb|ub] */
void load_crhmc_problem(SpMat &A, VT &b, VT &lb, VT &ub, int &dimension,
                        std::string problem_name) {
   {
    std::string fileName("./data/");
    fileName.append(problem_name);
    fileName.append(".mm");
    if(!exists_check(filename)){
      std::cerr<<"Problem does not exist.\n";
      exit(1);}
    SpMat X;
    loadMarket(X, fileName);
    int m = X.rows();
    dimension = X.cols() - 1;
    A = X.leftCols(dimension);
    b = VT(X.col(dimension));
  }
  {
    std::string fileName("./data/");
    fileName.append(problem_name);
    fileName.append("_bounds.mm");
    if(!exists_check(filename)){
      std::cerr<<"Problem does not exist.\n";
      exit(1);}
    SpMat bounds;
    loadMarket(bounds, fileName);
    lb = VT(bounds.col(0));
    ub = VT(bounds.col(1));
  }
}
template <int simdLen>
void run_main(std::string problem_name, int n_samples = 80000,
              int n_burns = 20000) {
  std::cerr<<"CRHMC on "<<problem_name<<"\n";
  using Solver =
      ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Input::Grad, simdLen>;
  RNG rng(1);
  Opts options;
  options.simdLen = simdLen;
  int dimension;
  SpMat A;
  VT b, lb, ub;
  load_crhmc_problem(A, b, lb, ub, dimension, problem_name);
  func_params params = func_params(Point(dimension), 0.5, 1);
  Func f(params);
  Grad g(params);
  Hess h(params);
  Input input = Input(dimension, f, g, h);
  input.Aeq = A;
  input.beq = b;
  input.lb = lb;
  input.ub = ub;
  std::ofstream stream;
  stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_" + problem_name +
              ".txt");
  std::cerr << "Finished loading data\n";
  options.EnableReordering = true;
  CrhmcProblem P = CrhmcProblem(input, options);
  stream<<"nnz = " <<P.Asp.nonZeros()<<"\n";
  std::cerr<<"nnz = " <<P.Asp.nonZeros()<<"\n";
  P.Asp.prune(1e-16);
  stream<<"nnz = " <<P.Asp.nonZeros()<<"\n";
  std::cerr<<"nnz = " <<P.Asp.nonZeros()<<"\n";
  std::cerr << "Finished Preparation process\n";
  P.print_preparation_time(stream);
  MT samples;
  sample<MT, CrhmcProblem, RNG, CRHMCWalk, NT, Point, Input, Solver, Opts>(
      samples, P, rng, n_samples, n_burns, input, options, stream);
  std::ofstream diagnostics_stream;
  diagnostics_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_" +
                          problem_name + "_diagnostics.txt");
  diagnose(samples, diagnostics_stream);
  std::ofstream samples_stream;
  samples_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_" +
                      problem_name + "_samples.txt");
  samples_stream << samples.transpose() << std::endl;
}
int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cerr
        << "Example Usage: ./crhmc_sample_sparse [problem_name] [simdLen] [n_samples] [n_burns]\n";
    std::cerr << "i.e.: ./crhmc_sample_sparse degen2 4 1000 500\n";
    exit(1);
  }
  if (atoi(argv[2]) == 1) {
    run_main<1>(argv[1], atoi(argv[3]), atoi(argv[4]));
  } else if (atoi(argv[2]) == 4) {
    run_main<4>(argv[1], atoi(argv[3]), atoi(argv[4]));
  } else if (atoi(argv[2]) == 8) {
    run_main<8>(argv[1], atoi(argv[3]), atoi(argv[4]));
  } else if (atoi(argv[2]) == 16) {
    run_main<16>(argv[1], atoi(argv[3]), atoi(argv[4]));
  }
  return 0;
}
