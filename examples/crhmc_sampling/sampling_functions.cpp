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
#include "volume/sampling_policies.hpp"
#include "ode_solvers/ode_solvers.hpp"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "sampling/random_point_generators.hpp"
#include "sampling/sampling.hpp"
#include "misc/misc.h"
#include "random.hpp"
#include <vector>
#include "random_walks/random_walks.hpp"
#include "generators/known_polytope_generators.h"


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
template <typename MT, typename VT, typename NT, typename StreamType>
void diagnose(MT &samples, StreamType &stream) {
  unsigned int min_ess = 0;
  print_diagnostics<NT, VT, MT>(samples, min_ess, stream);
  NT max_psrf = max_interval_psrf<NT, VT, MT>(samples);
  stream << "max_psrf: " << max_psrf << std::endl;
  stream << "min ess " << min_ess << std::endl;
}
using NT = double;
using Kernel = Cartesian<NT>;
using Point = typename Kernel::Point;
using Func = GaussianFunctor::FunctionFunctor<Point>;
using Grad = GaussianFunctor::GradientFunctor<Point>;
using Hess = GaussianFunctor::HessianFunctor<Point>;
using PolytopeType = HPolytope<Point>;
using MT = PolytopeType::MT;
using Input = crhmc_input<MT, Point, Func, Grad, Hess>;
using CrhmcProblem = crhmc_problem<Point, Input>;
using func_params = GaussianFunctor::parameters<NT, Point>;
using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;
template <int simdLen>
void run_main(std::string problem_name, int n_samples = 80000,
              int n_burns = 20000) {
  std::cerr << "CRHMC on " << problem_name << "\n";
  using Solver =
      ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Input::Grad, simdLen>;
  RNG rng(1);
  
  PolytopeType HP=generate_simplex<PolytopeType>(2,false);
  int dimension = HP.dimension();
  func_params params = func_params(Point(dimension), 0.5, 1);
  Func f(params);
  Grad g(params);
  Hess h(params);
  std::list<Point> PointList;
  crhmc_sampling<std::list<Point>, PolytopeType, RNG, CRHMCWalk, NT, Point, Grad, Func, Hess, Solver>(
      PointList, HP, rng, 1, n_samples, n_burns, g, f, h, simdLen);
  MT samples = MT(dimension, PointList.size());
  int i=0;
  for (std::list<Point>::iterator it = PointList.begin(); it != PointList.end(); ++it){
    samples.col(i) = (*it).getCoefficients();
    i++;
  }
  std::ofstream diagnostics_stream;
  diagnostics_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_" +
                          problem_name + "_diagnostics.txt");
  diagnose<MT, VT, NT, std::ofstream>(samples, diagnostics_stream);
  std::ofstream samples_stream;
  samples_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_" +
                      problem_name + "_samples.txt");
  samples_stream << samples.transpose() << std::endl;
}
int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cerr << "Example Usage: ./crhmc_sample_sparse [problem_name] "
                 "[simdLen] [n_samples] [n_burns]\n";
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
