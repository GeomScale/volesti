// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// VolEsti example

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "misc/misc.h"
#include "preprocess/crhmc/crhmcProblem.h"
#include "preprocess/crhmc/crhmc_input.h"
#include <fstream>
#include <time.h> /* clock_t, clock, CLOCKS_PER_SEC */
using NT=double;
using Kernel=Cartesian<NT>;
using Point=typename Kernel::Point;
using Hpolytope=HPolytope<Point>;
using CrhmcProblem = crhmcProblem<Point>;
using Input = crhmc_input<MT, NT>;
using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

void benchmark(std::string fileName) {
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
  double tstart;
  std::cout << "CRHMC polytope preparation for " << fileName << std::endl;
  tstart = (double)clock() / (double)CLOCKS_PER_SEC;
  CrhmcProblem P = CrhmcProblem(input);
  std::cout << "Preparation completed in time, ";
  std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart << " secs "
            << std::endl;
  std::cout << "The resulting matrix has " << P.Asp.nonZeros() << " nonZeros"
            << std::endl
            << std::endl;
}
int main() {

  std::cout << "CRHMC polytope preparation" << std::endl << std::endl;

  // Hpolytope HP = generate_cube<Hpolytope>(1000, false);
  int d = 100000;
  Input input = Input(d);
  input.lb = -VT::Ones(d);
  input.lb = VT::Ones(d);
  double tstart;

  std::cout << "CRHMC polytope preparation (cube-100000)" << std::endl;

  tstart = (double)clock() / (double)CLOCKS_PER_SEC;
  CrhmcProblem P = CrhmcProblem(input);

  std::cout << "Preparation completed in time, ";
  std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart << " secs "
            << std::endl;
  std::cout << "The resulting matrix has " << P.Asp.nonZeros() << " nonZeros"
            << std::endl
            << std::endl;
  benchmark("../metabolic_full_dim/polytope_e_coli.ine");
  benchmark("../netlib/afiro.ine");
  benchmark("../metabolic_full_dim/polytope_iAB_RBC_283.ine");
  benchmark("../metabolic_full_dim/polytope_recon1.ine");
  // benchmark("/content/drive/MyDrive/Polytopes/polytope_iAB_RBC_283.ine");

  return 0;
}
