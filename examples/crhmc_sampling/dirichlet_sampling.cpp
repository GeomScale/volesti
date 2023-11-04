// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2023 Vissarion Fisikopoulos
// Copyright (c) 2018-2023 Apostolos Chalkis
// Copyright (c) 2020-2023 Elias Tsigaridas

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with
// Riemannian Hamiltonian
// Monte Carlo in a Constrained Space"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
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
#include "helper_functions.hpp"
using NT = double;
using Kernel = Cartesian<NT>;
using Point = typename Kernel::Point;
using Func = DirichletFunctor::FunctionFunctor<Point>;
using Grad = DirichletFunctor::GradientFunctor<Point>;
using Hess = GaussianFunctor::HessianFunctor<Point>;
using PolytopeType = HPolytope<Point>;
using MT = PolytopeType::MT;
using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
using func_params = DirichletFunctor::parameters<NT, Point>;
using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;

template <int simdLen>
void sample_dirichlet(int n_samples = 80000,
              int n_burns = 20000){
  using SpMat = Eigen::SparseMatrix<NT>;
  using ConstraintProblem =constraint_problem<SpMat, Point>;
  using Input = crhmc_input<SpMat, Point, Func, Grad, Hess>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Solver =
      ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Input::Grad, simdLen>;
  
  int dim = 3;
  RNG rng(dim);
  MT A(1,dim);
  A <<  1.0, 1.0, 1.0;
  SpMat Aeq(1,3);
  Aeq.coeffRef(0,0) = NT(1);
  Aeq.coeffRef(0,1) = NT(1);
  Aeq.coeffRef(0,2) = NT(1);
  VT beq = VT::Ones(1), lb = VT::Zero(dim), ub = VT::Ones(dim);

  VT a_vec(dim);
  //a_vec << 2.0, 3.0, 4.0;
  a_vec << 2.0, 3.0, 4.0;
  
  std::string problem_name("dirichlet");
  ConstraintProblem problem = ConstraintProblem(dim);
  problem.set_equality_constraints(Aeq, beq);
  problem.set_bounds(lb, ub);
  func_params params = func_params(a_vec);
  Func *f = new Func(params);
  Grad *g = new Grad(params);
  Hess *h = NULL;
  std::list<Point> PointList;
  //crhmc_sampling<std::list<Point>, ConstraintProblem, RNG, CRHMCWalk, NT, Point, Grad, Func, Hess, Solver>(
  //    PointList, problem, rng, 1, n_samples, n_burns, g, f, h, simdLen);
  execute_crhmc<ConstraintProblem, RNG, std::list<Point>, Grad, Func, Hess, CRHMCWalk, simdLen>
            (problem, rng, PointList, 1, n_samples, n_burns, g, f, h);
  MT samples = MT(dim, PointList.size());
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

int main() {
  int n_samples = 5000;
  int n_burns = 2000;
            
  std::cerr<<"Sampling Dirichlet\n";
  sample_dirichlet<1>(n_samples, n_burns);
  sample_dirichlet<4>(n_samples, n_burns);
  sample_dirichlet<8>(n_samples, n_burns);
  sample_dirichlet<16>(n_samples, n_burns);

  return 0;
}
