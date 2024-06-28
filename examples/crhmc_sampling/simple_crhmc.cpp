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
#include "diagnostics/effective_sample_size.hpp"
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

template 
<
typename Polytope,
typename Point,
typename RandomNumberGenerator
>
void sample_aBW (Polytope &P, RandomNumberGenerator &rng, std::list<Point> &randPoints, unsigned int const&N)
{
        Point p = P.ComputeInnerBall().first;
        typedef typename AcceleratedBilliardWalk::template Walk
                <
                        Polytope,
                        RandomNumberGenerator
                > walk;

        typedef RandomPointGenerator <walk> RandomPointGenerator;
        PushBackWalkPolicy push_back_policy;

        RandomPointGenerator::apply(P, p, N, 1, randPoints,
                                    push_back_policy, rng);
}

double duration, maxPsrf;
unsigned int minEss;

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
    HP = random_orderpoly<PolytopeType, NT>(dim, m, 100);
    std::cout << "Sampling from Order Polytope" << std::endl;
  } else {
    HP = skinny_random_hpoly<PolytopeType, NT, PolyRNGType>(dim, m, true, NT(10000), 100);
    std::cout << "Sampling from Random skinny Polytope" << std::endl;
    // HP = generate_skinny_cube<PolytopeType>(20);
    // rotating<MT>(HP);
  }

  //HP.print();
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
    std::cout << "Using aBW walk" << std::endl;
    sample_aBW(HP, rng, PointList, n_samples);
  }

  stop = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> total_time = stop - start;
  std::cout << "Done in " << total_time.count() << '\n';
  duration = total_time.count();

  MT samples = MT(dim, PointList.size());
  int i=0;
  for (std::list<Point>::iterator it = PointList.begin(); it != PointList.end(); ++it){
    samples.col(i) = (*it).getCoefficients();
    i++;
  }
  unsigned int min_ess;
  NT max_psrf;
  effective_sample_size<NT, VT>(samples, min_ess);
  max_psrf = max_interval_psrf<NT,VT,MT>(samples);
  std::cerr << "min_ess: " << min_ess << '\n';
  std::cerr<<"max_psrf: "<< max_psrf <<"\n";
  minEss = min_ess;
  maxPsrf = max_psrf;
  /*
  std::ofstream samples_stream;
  samples_stream.open("CRHMC_SIMD_" + std::to_string(simdLen) + "_simplex" + "_samples.txt");
  samples_stream << samples.transpose() << std::endl;
  */
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
    if(argc == 9)
    {

      std::ofstream samples_stream;
      samples_stream.open("aoutput_yay" + std::to_string(atoi(argv[1])) + ".txt");
      for(int n = 20; n <= 100; n += 20)
      {
        std::cout << "\n\n" << n << '\n';
        int m = n * 10;
        double dur[2][2];
        double psrf[2][2];
        unsigned int ess[2][2];
        for(int a = 0; a <= 1; ++a)
          for(int b = 0; b <= 1; ++b) {
            if(atoi(argv[1]) == 8)
              run_main<8>(50*n, 0, n, m, a, b);
            else if(atoi(argv[1]) == 16)
              run_main<16>(50*n, 0, n, m, a, b);
            dur[a][b] = duration;
            psrf[a][b] = maxPsrf;
            ess[a][b] = minEss;
          }
        samples_stream << "\n\n";
        samples_stream << n << ' ' << m << "           order polytope          skinny polytope\n";
        samples_stream << "CRHMC: time:         " << dur[1][1] << "             " << dur[0][1] << '\n';
        samples_stream << "       psrf:         " << psrf[1][1] << "             " << psrf[0][1] << '\n';
        samples_stream << "       ess:          " << ess[1][1] << "                 " << ess[0][1] << '\n';
        samples_stream << '\n';
        samples_stream << "ABW:   time:         " << dur[1][0] << "             " << dur[0][0] << '\n';
        samples_stream << "       psrf:         " << psrf[1][0] << "             " << psrf[0][0] << '\n';
        samples_stream << "       ess:          " << ess[1][0] << "                 " << ess[0][0] << '\n';
        samples_stream << std::endl;
      }
      exit(1);
    }
    std::cerr << "Example Usage: ./simple_crhmc "
                 "[simdLen] [n_samples] [n_burns] [dimension] [facets] [if_order_poly] [if_crhmc_walk]\n";
    std::cerr << "i.e.: ./simple_crhmc 4 1000 500 20 150 1 0\n";
    exit(1);
  }
  // std::cerr << "To plot: python3 ../python_utilities/plot_samples.py <CRHMC_SIMD_4_simplex_samples.txt --save"<<"\n";

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
