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

template <typename NT,int simdLen>
void test_simdLen_sampling(int n_samples = 100000, int n_burns = -1,int dim=2,int walk_length=1,int burn_steps=1){
  std::cerr<<"--------------------------simdLen= "<<simdLen<<"\n";
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using pts = std::vector<Point>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;
  using Func = GaussianFunctor::FunctionFunctor<Point>;
  using Grad = GaussianFunctor::GradientFunctor<Point>;
  using Hess = GaussianFunctor::HessianFunctor<Point>;
  using func_params=GaussianFunctor::parameters<NT, Point>;
  using Input = crhmc_input<MT, Point, Func, Grad, Hess>;
  //using Input = crhmc_input<MT, Point>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Solver = ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad,simdLen>;
  using Opts = opts<NT>;
  RandomNumberGenerator rng(1);
  if (n_burns == -1) {
    n_burns = n_samples / 2;
  }
  func_params params = func_params(Point(dim), 4, 1);
  Func f(params);
  Grad g(params);
  Hess h(params);
  Opts options;
  options.simdLen=simdLen;
  //Input input =Input(dim);
  Input input = Input(dim, f, g, h);
  MT A = MT::Ones(5, dim);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = 10 * VT::Ones(5, 1);
  input.Aineq = A;
  input.bineq = b;
  CrhmcProblem P = CrhmcProblem(input, options);
  P.print();
  Point x0=Point(P.center);
  CRHMCWalk::parameters<NT, Grad> crhmc_params(g, dim, options);
  crhmc_params.eta = 0.2;
  crhmc_params.momentum = 0.8;
  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator, Grad, Func,
                  Solver>
      crhmc(P, x0, g, f, crhmc_params);
  MT samples = MT(dim, (n_samples - n_burns)*simdLen);
  #ifdef TIME_KEEPING
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::system_clock::now();
  #endif
  #ifdef TIME_KEEPING
    std::chrono::time_point<std::chrono::high_resolution_clock> start_file,
        end_file;
    std::chrono::duration<double> total_time_file =
        std::chrono::duration<double>::zero();
  #endif
    int j = 0;
    for (int i = 0; i < n_samples; i++) {
      if (i % 1000 == 0) {
        std::cerr << i << " out of " << n_samples << "\n";
      }
      for (int k = 0; k < burn_steps; k++) {
        crhmc.apply(rng, walk_length, true);
      }
  #ifdef TIME_KEEPING
      start_file = std::chrono::system_clock::now();
  #endif
      if (i >= n_burns) {
        MT sample = crhmc.getPoints();
        samples(Eigen::all,Eigen::seq((i-n_burns)*simdLen,(i-n_burns)*simdLen+simdLen-1)) = sample;
        j++;
      }
  #ifdef TIME_KEEPING
      end_file = std::chrono::system_clock::now();
      total_time_file += end_file - start_file;
  #endif
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
  #ifdef TIME_KEEPING
    start_file = std::chrono::system_clock::now();
  #endif
    std::cerr << "Writing samples in a file \n";
    std::cout << samples.transpose() << std::endl;
  #ifdef TIME_KEEPING
    end_file = std::chrono::system_clock::now();
    total_time_file += end_file - start_file;
    std::cerr << "Time for writing the file: " << total_time_file.count() << "\n";
  #endif
}

template void test_simdLen_sampling<double,4>(int n_samples = 100000, int n_burns = -1,int dim=2,int walk_length=1,int burn_steps=1);
template void test_simdLen_sampling<double,1>(int n_samples = 100000, int n_burns = -1,int dim=2,int walk_length=1,int burn_steps=1);
int main(int argc, char *argv[]) {

  if(argc != 5){
    std::cerr<< "Usage: ./crhmc_sampling [simdLen= 1 or 4] [n_samples] [n_burns] [dimension] \n";
  }else{
  if (atoi(argv[1])==1){
    test_simdLen_sampling<double,1>(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
  }
  else if (atoi(argv[1])==4){
    test_simdLen_sampling<double,4>(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
  }
  }
  return 0;
}
