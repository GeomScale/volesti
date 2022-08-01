// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with Riemannian Hamiltonian
//Monte Carlo in a Constrained Space"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "misc/misc.h"
#include "preprocess/crhmc/crhmcProblem.h"
#include "preprocess/crhmc/crhmc_input.h"
#include <fstream>
#include <iostream>
#include <vector>

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> HPOLYTOPE;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef crhmcProblem<Point> CRHMC_PROBLEM;
typedef crhmc_input<MT, NT> INPUT;

int main(int argc, char *argv[]) {
  unsigned d = 2;
  MT A = MT::Ones(5, d);
  A << 1, 0, -0.25, -1, 2.5, 1, 0.4, -1, -0.9, 0.5;
  VT b = 10 * VT::Ones(5, 1);
  HPOLYTOPE HP1 = HPolytope<Point>(d, A, b);
  std::cout << "Polytope HP1: \n";
  HP1.print();
  std::cout << "\n";
  INPUT input = INPUT(d);
  input.Aineq = A;
  input.bineq = b;
  CRHMC_PROBLEM P1 = CRHMC_PROBLEM(input);
  P1.print();

  std::string fileName("../../test/metabolic_full_dim/polytope_e_coli.ine");
  std::cout << "Reading input from file..." << std::endl;
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  inp.open(fileName, std::ifstream::in);
  read_pointset(inp, Pin);

  HPOLYTOPE HP2(Pin);
  // std::cout << "Polytope HP2: \n";
  // HP2.print();
  //   std::cout << "\n";
  CRHMC_PROBLEM P2 = CRHMC_PROBLEM(HP2);

  P2.print("coli_crhmc_polytope.txt");

  return 0;
}
