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
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "misc/misc.h"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include <fstream>
#include <iostream>
#include <time.h> /* clock_t, clock, CLOCKS_PER_SEC */
#include <vector>

using NT = double;
using Kernel = Cartesian<NT>;
using Point = Kernel::Point;
using PolytopeType = HPolytope<Point>;
using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
using Input = crhmc_input<MT, Point>;
using CrhmcProblem = crhmc_problem<Point, Input>;
using Opts = opts<NT>;

int main(int argc, char *argv[]) {
  unsigned d = 2;
  MT A = MT::Ones(5, d);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = 10 * VT::Ones(5, 1);
  PolytopeType HP1 = HPolytope<Point>(d, A, b);
  std::cout << "Polytope HP1: \n";
  HP1.print();
  std::cout << "\n";
  Input input = Input(d);
  input.Aineq = A;
  input.bineq = b;
  CrhmcProblem P1 = CrhmcProblem(input);
  P1.print();

  std::string fileName("../../test/metabolic_full_dim/polytope_e_coli.ine");
  std::cout << "Reading input from file..." << std::endl;
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  inp.open(fileName, std::ifstream::in);
  read_pointset(inp, Pin);

  PolytopeType HP2(Pin);
  // std::cout << "Polytope HP2: \n";
  // HP2.print();
  //   std::cout << "\n";
  d = HP2.dimension();
  Input input2 = Input(d);
  input2.Aineq = HP2.get_mat();
  input2.bineq = HP2.get_vec();
  Opts options;
  options.EnableReordering = false;
  double tstart = (double)clock() / (double)CLOCKS_PER_SEC;
  CrhmcProblem P2 = CrhmcProblem(input2, options);
  std::cout << "Preparation completed in time, ";
  std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart << std::endl;
  std::cout << "Number of nonZeros= " << P2.Asp.nonZeros() << "\n";
  // P2.print("coli_crhmc_polytope.txt");

  return 0;
}
