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
  if (argc != 2) {
    std::cout << "Usage: ./crhmc_prepare [file_name_string]\n";
    std::cout << "Example file name= "
                 "../../test/metabolic_full_dim/polytope_e_coli.ine\n";
    exit(1);
  }
  std::string fileName(argv[1]);
  std::cout << "Reading input from file..." << std::endl;
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  inp.open(fileName, std::ifstream::in);
  read_pointset(inp, Pin);

  PolytopeType HP(Pin);
  int d = HP.dimension();
  Input input = Input(d);
  input.Aineq = HP.get_mat();
  input.bineq = HP.get_vec();
  Opts options;
  options.EnableReordering = false;
  double tstart = (double)clock() / (double)CLOCKS_PER_SEC;
  CrhmcProblem P = CrhmcProblem(input, options);
  std::cout << "Preparation completed in time, ";
  std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart << std::endl;
  std::cout << "Number of nonZeros= " << P.Asp.nonZeros() << "\n";
  return 0;
}
