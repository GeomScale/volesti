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
#include <assert.h>
#include <fstream>
#include <chrono>
using NT = double;
using Kernel = Cartesian<NT>;
using Point = typename Kernel::Point;
using Hpolytope = HPolytope<Point>;
using Input = crhmc_input<MT, Point>;
using CrhmcProblem = crhmc_problem<Point, Input>;
using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

inline bool exists_check(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

double benchmark(std::string fileName) {
  std::ifstream inp;
  inp.open(fileName, std::ifstream::in);
  std::vector<std::vector<NT>> Pin;
  read_pointset(inp, Pin);
  inp.close();
  Hpolytope HP(Pin);
  int d = HP.dimension();
  Input input = Input(d);
  input.Aineq = HP.get_mat();
  input.bineq = HP.get_vec();
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::cout << "CRHMC polytope preparation for " << fileName << std::endl;

  start = std::chrono::system_clock::now();
  CrhmcProblem P = CrhmcProblem(input);
  end = std::chrono::system_clock::now();

  std::cout << "Preparation completed in time, ";
  std::chrono::duration<double> elapsed_seconds = end - start;
  double preparation_time = elapsed_seconds.count();
  std::cout << preparation_time << " secs " << std::endl;
  std::cout << "The resulting matrix has " << P.Asp.nonZeros() << " nonZeros"
            << std::endl;
  return preparation_time;
}

int main() {

  std::cout
      << "---------------CRHMC polytope preparation benchmarking---------------"
      << std::endl
      << std::endl;

  int d = 100000;
  Input input = Input(d);
  input.lb = -VT::Ones(d);
  input.ub = VT::Ones(d);
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

  std::cout << "CRHMC polytope preparation 100000 dimensional Cube "
            << std::endl;

  start = std::chrono::system_clock::now();
  CrhmcProblem P = CrhmcProblem(input);
  end = std::chrono::system_clock::now();

  std::cout << "Preparation completed in time, ";
  std::chrono::duration<double> elapsed_seconds = end - start;
  double preparation_time = elapsed_seconds.count();
  std::cout << preparation_time << " secs " << std::endl;
  std::cout << "The resulting matrix has " << P.Asp.nonZeros() << " nonZeros"
            << std::endl;
  assert(preparation_time < 0.5);
  std::cout << "Assertion (preparation_time< 0.5 secs) passed!" << std::endl
            << std::endl;
  if (exists_check("../test/metabolic_full_dim/polytope_e_coli.ine")) {
    preparation_time =
        benchmark("../test/metabolic_full_dim/polytope_e_coli.ine");
    assert(preparation_time < 2.0);
    std::cout << "Assertion (preparation_time< 2 secs) passed!" << std::endl
              << std::endl;
  }
  if (exists_check("../test/netlib/afiro.ine")) {
    preparation_time = benchmark("../test/netlib/afiro.ine");
    assert(preparation_time < 0.3);
    std::cout << "Assertion (preparation_time< 0.3 secs) passed!" << std::endl
              << std::endl;
  }
  if (exists_check("../test/metabolic_full_dim/polytope_iAB_RBC_283.ine")) {
    preparation_time =
        benchmark("../test/metabolic_full_dim/polytope_iAB_RBC_283.ine");
    assert(preparation_time < 400);
    std::cout << "Assertion (preparation_time< 400 secs) passed!" << std::endl
              << std::endl;
  }
  return 0;
}
