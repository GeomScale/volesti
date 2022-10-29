// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <pthread.h>
#include <string>
#include <typeinfo>
#include <unistd.h>
#include <vector>

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "doctest.h"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "ode_solvers/oracle_functors.hpp"
#include "generators/known_polytope_generators.h"
#include "convex_bodies/hpolytope.h"
#include "misc/misc.h"

template <typename NT> void test_crhmc_polytope_preprocessing() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using Input = crhmc_input<MT, Point>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using PolytopeType = HPolytope<Point>;
  using Opts = opts<NT>;

  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  std::string fileName("../test/metabolic_full_dim/polytope_e_coli.ine");
  inp.open(fileName, std::ifstream::in);
  read_pointset(inp, Pin);
  inp.close();
  PolytopeType HP(Pin);
  Opts options = Opts();
  int d = HP.dimension();
  Input input = Input(d);
  input.Aineq = HP.get_mat();
  input.bineq = HP.get_vec();
  options.EnableReordering = false;
  CrhmcProblem P = CrhmcProblem(input, options);

  int m = 342;
  int n = 366;
  std::ifstream testdata;
  std::string testDataFileName("../test/crhmc_polytope_test_output.txt");
  testdata.open(testDataFileName, std::ifstream::in);
  int size;
  testdata >> size;

  // CHECK(size == m);
  testdata >> size;

  // CHECK(size == n);
  MT A = MT(P.Asp);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      NT Matrxidata;
      testdata >> Matrxidata;
      CHECK(std::abs(A(i, j) - (Matrxidata)) < 0.001);
    }
  }
}
template <typename NT> void test_crhmc_fixed_var_polytope() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using Input = crhmc_input<MT, Point>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using PolytopeType = HPolytope<Point>;
  using Opts = opts<NT>;
  unsigned d = 2;
  MT A = MT(5, d);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = VT(5, 1);
  b << 10, 10, 10, 10, 10;
  Input input = Input(d);
  input.Aineq = A;
  input.bineq = b;
  input.lb(1) = 10;
  input.ub(1) = 10;
  Opts options = Opts();
  options.EnableReordering = false;
  CrhmcProblem P = CrhmcProblem(input, options);
  MT Aout = MT(P.Asp);
  MT Acheck = MT(5, 6);
  VT a = VT(5);
  a << 0.5, -0.25, 0.625, 0.2, -0.45;
  Acheck << a, MT::Identity(5, 5);
  int n = Aout.cols();
  int m = Aout.rows();
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      CHECK(std::abs(Aout(i, j) - Acheck(i, j)) < 0.001);
    }
  }
}
template <typename NT> void test_crhmc_dependent_polytope() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using Input = crhmc_input<MT, Point>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using SpMat = typename CrhmcProblem::SpMat;
  using InputSparse = crhmc_input<SpMat, Point>;
  using PolytopeType = HPolytope<Point>;
  using Opts = opts<NT>;
  unsigned d = 3;
  MT Aeq = MT(3, d);
  Aeq << 1, 2, 3, 1, 1, 1, 2, 3, 4;
  VT beq = VT(3, 1);
  beq << 10, 10, 20;
  Input input = Input(d);
  input.Aeq = Aeq;
  input.beq = beq;
  Opts options = Opts();
  options.EnableReordering = true;
  CrhmcProblem P = CrhmcProblem(input, options);
  CHECK(P.equations() == 2);

  SpMat A = Aeq.sparseView();
  InputSparse input_sparse = InputSparse(d);
  input_sparse.Aeq = A;
  input_sparse.beq = beq;
  crhmc_problem<Point, InputSparse> Q = crhmc_problem<Point, InputSparse>(input_sparse,options);
  CHECK(Q.equations() == 2);
}

template <typename NT> void test_center_computation() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using Func = GaussianFunctor::FunctionFunctor<Point>;
  using Grad = GaussianFunctor::GradientFunctor<Point>;
  using Hess = GaussianFunctor::HessianFunctor<Point>;
  using func_params=GaussianFunctor::parameters<NT, Point>;
  using Input = crhmc_input<MT, Point, Func, Grad, Hess>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using PolytopeType = HPolytope<Point>;
  using Opts = opts<NT>;
  using PolytopeType = HPolytope<Point>;
  unsigned dim = 2;
  Point mean=Point(VT::Ones(dim));
  mean=mean*0.5;
  func_params params = func_params(mean, 0.5, 1);
  Func f(params);
  Grad g(params);
  Hess h(params);
  Opts options = Opts();
  options.EnableReordering = true;
  PolytopeType HP = generate_cross<PolytopeType>(2, false);
  Input input = Input(dim, f, g, h);
  input.Aineq = HP.get_mat();
  input.bineq = HP.get_vec();
  CrhmcProblem P = CrhmcProblem(input, options);
  VT analytic_ctr=VT(P.dimension());
  analytic_ctr<<  0.0970, 0.0970, 1.1939, 1.0000, 1.0000, 0.8061;
  VT lewis_center=VT(P.dimension());
  lewis_center<< -0.0585, -0.0585, -0.1171, 0, 0, 0.1171;
  for(int i=0;i<dim;i++){
    CHECK(std::abs(P.analytic_ctr(i) - analytic_ctr(i)) < 0.001);
    CHECK(std::abs(P.center(i) - (lewis_center(i))) < 0.001);
  }

}

template <typename NT> void test_unbounded() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using Input = crhmc_input<MT, Point>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Opts = opts<NT>;
  Input input= Input(1);
  Opts options = Opts();
  CrhmcProblem P = CrhmcProblem(input,options);
  CHECK(P.terminate);

}
template <typename NT> void call_test_crhmc_fixed_var_polytope(){
  std::cout << "--- Testing fixed vars" << std::endl;
  test_crhmc_fixed_var_polytope<NT>();
}
template <typename NT>void call_test_crhmc_dependent_polytope(){
  std::cout << "--- Testing dep vars" << std::endl;
  test_crhmc_dependent_polytope<NT>();
}
template <typename NT> void call_test_crhmc_preprocesssing() {
  std::cout << "--- Testing CRHMC data preprocessing" << std::endl;
  test_crhmc_polytope_preprocessing<NT>();
}
template <typename NT> void call_test_center_computation(){
  std::cout << "--- Testing CRHMC polytope-center computation" << std::endl;
  test_center_computation<NT>();
}

template <typename NT> void call_test_invalid_polytopes(){
  std::cout << "--- Testing CRHMC on Unbounded" << std::endl;
  test_unbounded<NT>();
}

TEST_CASE("test_preparation_crhmc") {
  call_test_crhmc_preprocesssing<double>();
}

TEST_CASE("test_fixed_vars_crhmc") {
  call_test_crhmc_fixed_var_polytope<double>();
}

TEST_CASE("test_dep_vars_crhmc") {
  call_test_crhmc_dependent_polytope<double>();
}

TEST_CASE("test_center_computation"){
  call_test_center_computation<double>();
}

TEST_CASE("test_invalid_polytopes"){
  call_test_invalid_polytopes<double>();
}
