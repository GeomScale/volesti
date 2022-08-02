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
#include "preprocess/crhmc/crhmcProblem.h"
#include "preprocess/crhmc/crhmc_input.h"

#include "convex_bodies/hpolytope.h"
#include "misc/misc.h"

template <typename NT> void test_crhmc_polytope_preprocessing() {
  typedef Cartesian<NT> Kernel;
  typedef typename Kernel::Point Point;
  typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
  typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
  typedef crhmcProblem<Point> CRHMC_PROBLEM;
  typedef crhmc_input<MT, NT> INPUT;
  typedef HPolytope<Point> HPOLYTOPE;

  /*
  unsigned d = 2;
  unsigned m=5;
  MT A = MT::Ones(m, d);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = 10 * VT::Ones(m, 1);
  INPUT input = INPUT(d);
  input.Aineq = A;
  input.bineq = b;
  CRHMC_PROBLEM P = CRHMC_PROBLEM(input);

  CHECK(std::abs(P.y(0) - (-4.8342)) < 0.001);
  CHECK(std::abs(P.center(0) - 3.75986) < 0.001);
*/
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  std::string fileName("../metabolic_full_dim/polytope_e_coli.ine");
  inp.open(fileName, std::ifstream::in);
  read_pointset(inp, Pin);

  HPOLYTOPE HP(Pin);
  d=HP2.dimension();
  INPUT input = INPUT(d);
  input.Aineq=HP2.get_mat();
  input.bineq=HP2.get_vec();
  CRHMC_PROBLEM P = CRHMC_PROBLEM(input);
  //CRHMC_PROBLEM P = CRHMC_PROBLEM(HP);

  int m = 174;
  int n = 198;
  std::ifstream testdata;
  std::string testDataFileName("../../examples/crhmc_prepare/outputMatrix.txt");
  testdata.open(testDataFileName, std::ifstream::in);
  int size;
  testdata >> size;

  CHECK(size == m);
  testdata >> size;

  CHECK(size == n);
  MT A=MT(P.Asp);
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
