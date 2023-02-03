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
#include "common.hpp"

template <int simdLen>
void sample_hpoly(int n_samples = 80000,
              int n_burns = 20000, int dim = 2) {
  using NT = double;
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using Func = ZeroScalarFunctor<Point>;
  using Grad = ZeroFunctor<Point>;
  using Hess = ZeroFunctor<Point>;
  using PolytopeType = HPolytope<Point>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using MT = PolytopeType::MT;
  using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;
  std::string problem_name("simplex");
  std::cerr << "CRHMC on " << problem_name << "\n";
  RNG rng(1);
  PolytopeType HP=generate_simplex<PolytopeType>(dim,false);
  int dimension = HP.dimension();
  Func * f = new Func;
  Grad * g = new Grad;
  std::list<Point> PointList;
  execute_crhmc< PolytopeType, RNG, std::list<Point>, Grad, Func, Hess, CRHMCWalk, simdLen>(
      HP, rng, PointList, 1, n_samples, n_burns, g, f);
  MT samples = MT(dimension, PointList.size());
  int i=0;
  for (std::list<Point>::iterator it = PointList.begin(); it != PointList.end(); ++it){
    samples.col(i) = (*it).getCoefficients();
    i++;
  }
  std::cerr<<"max_psrf: "<< max_interval_psrf<NT,VT,MT>(samples)<<"\n";
  std::ofstream samples_stream;
  samples_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_" +
                      problem_name + "_samples.txt");
  samples_stream << samples.transpose() << std::endl;
}

template<int simdLen>
void run_main(int n_samples = 80000,
              int n_burns = 20000,
              int dimension = 2){
  std::cerr<<"Sampling HPolytope\n";
  sample_hpoly<simdLen>(n_samples, n_burns, dimension);
}
int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cerr << "Example Usage: ./simple_crhmc "
                 "[simdLen] [n_samples] [n_burns] [dimension]\n";
    std::cerr << "i.e.: ./simple_crhmc 4 1000 500 2\n";
    exit(1);
  }
  std::cerr << "To plot: python3 ../python_utilities/plot_samples.py <CRHMC_SIMD_4_simplex_samples.txt --save"<<"\n";
  if (atoi(argv[1]) == 1) {
    run_main<1>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
  } else if (atoi(argv[1]) == 4) {
    run_main<4>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
  } else if (atoi(argv[1]) == 8) {
    run_main<8>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
  } else if (atoi(argv[1]) == 16) {
    run_main<16>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
  }
  return 0;
}
