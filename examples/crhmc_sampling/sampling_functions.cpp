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
#include <unsupported/Eigen/SparseExtra>


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
using func_params = GaussianFunctor::parameters<NT, Point>;
using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;
template <int simdLen>
void sample_hpoly(int n_samples = 80000,
              int n_burns = 20000) {
  std::string problem_name("simplex");
  std::cerr << "CRHMC on " << problem_name << "\n";
  using Input = crhmc_input<MT, Point, Func, Grad, Hess>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
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
inline bool exists_check(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}
/*Problem on the form X=[A|b] bounds=[lb|ub] */
void load_crhmc_problem(SpMat &A, VT &b, VT &lb, VT &ub, int &dimension,
                        std::string problem_name) {
   {
    std::string fileName("./data/");
    fileName.append(problem_name);
    fileName.append(".mm");
    if(!exists_check(fileName)){
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
    if(!exists_check(fileName)){
      std::cerr<<"Problem does not exist.\n";
      exit(1);}
    SpMat bounds;
    loadMarket(bounds, fileName);
    lb = VT(bounds.col(0));
    ub = VT(bounds.col(1));
  }
}
template <int simdLen>
void sample_sparse_problem(int n_samples = 80000,
              int n_burns = 20000){
  using SpMat = Eigen::SparseMatrix<NT>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using ConstraintProblem =constraint_problem<SpMat, Point>;
  std::string problem_name("afiro");
  std::cerr << "CRHMC on " << problem_name << "\n";
  using Input = crhmc_input<SpMat, Point, Func, Grad, Hess>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Solver =
      ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Input::Grad, simdLen>;

  RNG rng(1);
  SpMat A;
  VT b, lb, ub;
  int dimension;
  load_crhmc_problem(A, b, lb, ub, dimension, problem_name);
  ConstraintProblem problem = ConstraintProblem(dimension);
  problem.set_equality_constraints(A, b);
  problem.set_bounds(lb, ub);
  func_params params = func_params(Point(dimension), 0.5, 1);
  Func f(params);
  Grad g(params);
  Hess h(params);
  std::list<Point> PointList;
  crhmc_sampling<std::list<Point>, ConstraintProblem, RNG, CRHMCWalk, NT, Point, Grad, Func, Hess, Solver>(
      PointList, problem, rng, 1, n_samples, n_burns, g, f, h, simdLen);
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
template<int simdLen>
void run_main(int n_samples = 80000,
              int n_burns = 20000){
  std::cerr<<"Sampling HPolytope\n";
  sample_hpoly<simdLen>(n_samples, n_burns);
  std::cerr<<"Sampling Sparse Problem\n";
  sample_sparse_problem<simdLen>(n_samples, n_burns);
}
int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cerr << "Example Usage: ./crhmc_sample_sparse "
                 "[simdLen] [n_samples] [n_burns]\n";
    std::cerr << "i.e.: ./crhmc_sample_ 4 1000 500\n";
    exit(1);
  }
  if (atoi(argv[1]) == 1) {
    run_main<1>(atoi(argv[2]), atoi(argv[3]));
  } else if (atoi(argv[1]) == 4) {
    run_main<4>(atoi(argv[2]), atoi(argv[3]));
  } else if (atoi(argv[1]) == 8) {
    run_main<8>(atoi(argv[2]), atoi(argv[3]));
  } else if (atoi(argv[1]) == 16) {
    run_main<16>(atoi(argv[2]), atoi(argv[3]));
  }
  return 0;
}
