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
#include "generators/h_polytopes_generator.h"
#include "generators/order_polytope_generator.h"
#include "volume/sampling_policies.hpp"
#include "ode_solvers/ode_solvers.hpp"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "sampling/random_point_generators.hpp"
#include "sampling/sampling.hpp"
#include "misc/misc.h"
#include "misc/poset.h"
#include "random.hpp"
#include <vector>
#include "random_walks/random_walks.hpp"
#include "generators/known_polytope_generators.h"
#include "helper_functions.hpp"
#include "volume/rotating.hpp"
#include "convex_bodies/orderpolytope.h"


template 
<
typename Polytope,
typename Point,
typename RandomNumberGenerator
>
void sample_cdhr (Polytope &P, RandomNumberGenerator &rng, std::list<Point> &randPoints, unsigned int const&N)
{
        Point p = P.ComputeInnerBall().first;
        typedef typename CDHRWalk::template Walk
                <
                        Polytope,
                        RandomNumberGenerator
                > walk;

        typedef RandomPointGenerator <walk> RandomPointGenerator;
        PushBackWalkPolicy push_back_policy;

        RandomPointGenerator::apply(P, p, N, 1, randPoints,
                                    push_back_policy, rng);
}


template <int simdLen>
void sample_hpoly(int n_samples = 80000,
              int n_burns = 20000, int dim = 20, int m = 150, bool order_poly = false, bool crhmc_walk = false) {
  using NT = double;
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using Func = ZeroScalarFunctor<Point>;
  using Grad = ZeroFunctor<Point>;
  using Hess = ZeroFunctor<Point>;
  using PolytopeType = HPolytope<Point>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using MT = PolytopeType::MT;
  typedef boost::mt19937 PolyRNGType;
  using RNG = BoostRandomNumberGenerator<boost::mt19937, NT>;


  RNG rng(dim);
  PolytopeType HP;
  if(order_poly) {
    HP = random_orderpoly<PolytopeType, NT>(dim, m);
    std::cout << "Sampling from Order Polytope" << std::endl;
  } else {
    HP = skinny_random_hpoly<PolytopeType, NT, PolyRNGType>(dim, m, false, NT(4000));
    std::cout << "Sampling from Random skinny Polytope" << std::endl;
    // HP = generate_skinny_cube<PolytopeType>(20);
    // rotating<MT>(HP);
  }

  // HP.print();
  Func * f = new Func;
  Grad * g = new Grad;
  std::list<Point> PointList;

  
  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
  start = std::chrono::high_resolution_clock::now();

  if(crhmc_walk) {
    std::cout << "Using CRHMC walk" << std::endl;
    execute_crhmc< PolytopeType, RNG, std::list<Point>, Grad, Func, Hess, CRHMCWalk, simdLen>(
      HP, rng, PointList, 1, n_samples, n_burns, g, f);
  } else {
    std::cout << "Using CDHR walk" << std::endl;
    sample_cdhr(HP, rng, PointList, n_samples);
  }

  stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> total_time = stop - start;
  std::cout << "Done in " << total_time.count() << '\n';

  MT samples = MT(dim, PointList.size());
  int i=0;
  for (std::list<Point>::iterator it = PointList.begin(); it != PointList.end(); ++it){
    samples.col(i) = (*it).getCoefficients();
    i++;
  }
  std::cerr<<"max_psrf: "<< max_interval_psrf<NT,VT,MT>(samples)<<"\n";
  std::ofstream samples_stream;
  samples_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_simplex" + "_samples.txt");
  samples_stream << samples.transpose() << std::endl;
  delete f;
  delete g;
}


template<int simdLen>
void run_main(int n_samples = 80000,
              int n_burns = 20000,
              int dimension = 20, int m = 150, bool order_poly = false, bool crhmc_walk = false){
  sample_hpoly<simdLen>(n_samples, n_burns, dimension, m, order_poly, crhmc_walk);
}
int main(int argc, char *argv[]) {
  if (argc != 8) {
    std::cerr << "Example Usage: ./simple_crhmc "
                 "[simdLen] [n_samples] [n_burns] [dimension] [facets] [if_order_poly] [if_crhmc_walk]\n";
    std::cerr << "i.e.: ./simple_crhmc 4 1000 500 20 150 1 0\n";
    exit(1);
  }
  std::cerr << "To plot: python3 ../python_utilities/plot_samples.py <CRHMC_SIMD_4_simplex_samples.txt --save"<<"\n";
  if (atoi(argv[1]) == 1) {
    run_main<1>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
  } else if (atoi(argv[1]) == 4) {
    run_main<4>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
  } else if (atoi(argv[1]) == 8) {
    run_main<8>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
  } else if (atoi(argv[1]) == 16) {
    run_main<16>(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
  }
  return 0;
}
