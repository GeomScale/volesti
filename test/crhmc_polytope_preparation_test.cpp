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
#include "preprocess/crhmc/crhmc_problem.h"
#include "preprocess/crhmc/crhmc_input.h"

#include "convex_bodies/hpolytope.h"
#include "misc/misc.h"

template <typename NT> void test_crhmc_polytope_preprocessing() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using CrhmcProblem = crhmc_problem<Point>;
  using Input = crhmc_input<MT, NT>;
  using PolytopeType = HPolytope<Point>;
  using Opts = opts<NT>;

  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  std::string fileName("../test/metabolic_full_dim/polytope_e_coli.ine");
  inp.open(fileName, std::ifstream::in);
  read_pointset(inp, Pin);
  inp.close();
  PolytopeType HP(Pin);
  int d = HP.dimension();
  Input input = Input(d);
  input.Aineq = HP.get_mat();
  input.bineq = HP.get_vec();
  Opts options;
  options.EnableReordering = false;
  CrhmcProblem P = CrhmcProblem(input, options);

  int m = 342;
  int n = 366;
  std::ifstream testdata;
  std::string testDataFileName("../examples/crhmc_prepare/outputMatrix.txt");
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

template <typename NT> void call_test_crhmc_preprocesssing() {
  std::cout << "--- Testing CRHMC data preprocessing" << std::endl;
  test_crhmc_polytope_preprocessing<NT>();
}

TEST_CASE("simple_example_crhmc") { call_test_crhmc_preprocesssing<double>(); }
