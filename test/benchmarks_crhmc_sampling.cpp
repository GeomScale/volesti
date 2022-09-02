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
#include "diagnostics/diagnostics.hpp"
#include "generators/known_polytope_generators.h"
#include "misc/misc.h"
#include "ode_solvers.hpp"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "random.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_int.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include <assert.h>
#include <chrono>
#include <fstream>
#include "preprocess/svd_rounding.hpp"
template <typename Polytope, typename NT>
Polytope read_polytope(std::string filename) {
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  inp.open(filename, std::ifstream::in);
  read_pointset(inp, Pin);
  Polytope P(Pin);
  return P;
}
template <typename NT> struct SimulationStats {
  std::string method;
  unsigned int walk_length;
  unsigned int min_ess = 0;
  NT max_psrf = NT(0);
  NT time_per_draw = NT(0);
  NT time_per_independent_sample = NT(0);
  NT average_acceptance_prob = NT(0);
  NT step_size = NT(0);

  friend std::ostream &operator<<(std::ostream &out,
                                  const SimulationStats &stats) {
    out << stats.method << "," << stats.walk_length << "," << stats.min_ess
        << "," << stats.max_psrf << "," << stats.time_per_draw << ","
        << stats.time_per_independent_sample << ","
        << stats.average_acceptance_prob << ","
        << "," << stats.step_size << std::endl;
    return out;
  }
};
struct InnerBallFunctor {

  // Gaussian density centered at the inner ball center
  template <typename NT, typename Point> struct parameters {
    unsigned int order;
    NT L;     // Lipschitz constant for gradient
    NT m;     // Strong convexity constant
    NT kappa; // Condition number
    NT R0;
    NT sigma;
    Point x0;

    parameters(Point x0_, NT R0_)
        : order(2), L(1), m(1), kappa(1), x0(x0_), R0(R0_), sigma(1){};
  };

  template <typename Point> struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT, Point> &params;

    GradientFunctor(parameters<NT, Point> &params_) : params(params_){};

    // The index i represents the state vector index
    Point operator()(unsigned int const &i, pts const &xs, NT const &t) const {
      if (i == params.order - 1) {
        Point y = (-1.0 / pow(params.sigma, 2)) * (xs[0] - params.x0);
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }
    Point operator()(Point const &x) const {
      Point y = (-1.0 / pow(params.sigma, 2)) * (x - params.x0);
      return y;
    }
  };

  template <typename Point> struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_){};

    // The index i represents the state vector index
    NT operator()(Point const &x) const {
      Point y = x - params.x0;
      return 1.0 / (2 * pow(params.sigma, 2)) * y.dot(y);
    }
  };
  template <typename Point> struct Hess {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;
    Hess(parameters<NT, Point> &params_) : params(params_){};

    Point operator()(Point const &x) const {
      return (1.0 / pow(params.sigma, 2)) * Point::all_ones(x.dimension());
    }
  };
};

inline bool exists_check(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}
template <typename NT, typename VT, typename MT>
NT check_interval_psrf(MT &samples, NT target = NT(1.2)) {
  NT max_psrf = NT(0);
  VT intv_psrf = interval_psrf<VT, NT, MT>(samples);
  unsigned int d = intv_psrf.rows();
  for (unsigned int i = 0; i < d; i++) {
    assert(intv_psrf(i) < target);
    if (intv_psrf(i) > max_psrf)
      max_psrf = intv_psrf(i);
  }
  return max_psrf;
}
template <typename NT, typename Polytope>
std::vector<SimulationStats<NT>> benchmark_polytope_sampling(
    Polytope &P, NT eta = NT(-1), unsigned int walk_length = 1,
    double target_time = std::numeric_limits<NT>::max(), bool rounding = false,
    bool centered = false, unsigned int max_draws = 80000,
    unsigned int num_burns = 20000) {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using MT = typename Polytope::MT;
  using VT = typename Polytope::VT;
  using Func = InnerBallFunctor::FunctionFunctor<Point>;
  using Grad = InnerBallFunctor::GradientFunctor<Point>;
  using Hess = InnerBallFunctor::Hess<Point>;
  using Input = crhmc_input<MT, Point, Func, Grad, Hess>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Opts = opts<NT>;
  using Solver = ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad>;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;

  SimulationStats<NT> rdhr_stats;
  SimulationStats<NT> crhmc_stats;

  std::pair<Point, NT> inner_ball;
  if (centered) {
    inner_ball.first = Point(P.dimension());
    inner_ball.second = NT(1); // dummy radius (not correct one)
  } else {
    inner_ball = P.ComputeInnerBall();
  }

  // Random number generator
  RandomNumberGenerator rng(1);

  // Chebyshev center
  Point x0 = inner_ball.first;
  NT R0 = inner_ball.second;
  unsigned int dim = x0.dimension();

  if (rounding) {
    std::cout << "SVD Rounding" << std::endl;
    svd_rounding<AcceleratedBilliardWalk, MT, VT>(P, inner_ball, walk_length,
                                                  rng);
  }
  // Declare oracles
  InnerBallFunctor::parameters<NT, Point> params(x0, R0);

  Grad F(params);
  Func f(params);
  Hess H(params);
  int max_actual_draws = max_draws - num_burns;
  unsigned int min_ess = 0;

  MT samples;
  samples.resize(dim, max_actual_draws);
  NT ETA;
  NT max_psrf;

  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
  Opts options;
  CRHMCWalk::parameters<NT, Grad> crhmc_params(F, dim, options);
  Input input = Input(P.dimension(), f, F, H);
  input.Aineq = P.get_mat();
  input.bineq = P.get_vec();

  CrhmcProblem crhmc_problem = CrhmcProblem(input);
  Point x_start(crhmc_problem.center);
  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator, Grad, Func,
                  Solver>
      crhmc(crhmc_problem, x_start, F, f, crhmc_params);

  min_ess = 0;

  std::cout
      << "Constrained Riemannian Hamiltonian Monte Carlo (Gaussian Density)"
      << std::endl;

  if (eta > 0)
    crhmc.solver->eta = eta;

  std::cout << "Burn-in" << std::endl;

  for (unsigned int i = 0; i < num_burns; i++) {
    if (i % 1000 == 0)
      std::cout << ".";
    crhmc.apply(rng, 1);
  }

  std::cout << std::endl;
  std::cout << "Sampling" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for (unsigned int i = 0; i < max_actual_draws; i++) {
    for (int k = 0; k < walk_length; k++) {
      crhmc.apply(rng, 1);
    }
    samples.col(i) = crhmc.getPoint().getCoefficients();
    if (i % 1000 == 0 && i > 0)
      std::cout << ".";
  }
  stop = std::chrono::high_resolution_clock::now();

  ETA = (NT)std::chrono::duration_cast<std::chrono::microseconds>(stop - start)
            .count();

  std::cout << std::endl;
#ifdef TIME_KEEPING
  std::chrono::duration<double> total_time = stop - start;
  std::cerr << "Total time: " << total_time.count() << "\n";
  assert(total_time.count() < target_time);
  std::cout << "Assertion (preparation_time< " << target_time
            << " secs) passed!" << std::endl
            << std::endl;
  crhmc.print_timing_information();
#endif
  print_diagnostics<NT, VT, MT>(samples, min_ess, std::cout);
  std::cout << "min ess " << min_ess << "us" << std::endl;
  std::cout << "Average time per sample: " << ETA / max_actual_draws << "us"
            << std::endl;
  std::cout << "Average time per independent sample: " << ETA / min_ess << "us"
            << std::endl;
  std::cout << "Step size (final): " << crhmc.solver->eta << std::endl;
  std::cout << "Discard Ratio: " << crhmc.discard_ratio << std::endl;
  std::cout << "Average Acceptance Probability: "
            << crhmc.average_acceptance_prob << std::endl;
  max_psrf = check_interval_psrf<NT, VT, MT>(samples);
  std::cout << "max_psrf: " << max_psrf << std::endl;
  std::cout << std::endl;

  crhmc_stats.method = "CRHMC";
  crhmc_stats.walk_length = walk_length;
  crhmc_stats.min_ess = min_ess;
  crhmc_stats.max_psrf = max_psrf;
  crhmc_stats.time_per_draw = ETA / max_actual_draws;
  crhmc_stats.time_per_independent_sample = ETA / min_ess;
  crhmc_stats.step_size = crhmc.solver->eta;
  crhmc_stats.average_acceptance_prob = crhmc.average_acceptance_prob;

  return std::vector<SimulationStats<NT>>{rdhr_stats, crhmc_stats};
}

template <typename NT, typename Point, typename HPolytope>
void test_benchmark_polytope(
    HPolytope &P, std::string &name, bool centered,
    double target_time = std::numeric_limits<NT>::max(), int walk_length = 1) {
  std::cout << "CRHMC polytope preparation for " << name << std::endl;
  std::vector<SimulationStats<NT>> results;
  NT step_size = 0;
  std::pair<Point, NT> inner_ball;
  std::ofstream outfile;
  std::cout << name << std::endl;
  outfile.open("results_" + name + "_new.txt");
  P.normalize();
  inner_ball = P.ComputeInnerBall();
  step_size = inner_ball.second / 10;
  results = benchmark_polytope_sampling<NT, HPolytope>(P, step_size, walk_length, target_time,
                                        false, centered);
  outfile << results[0];
  outfile << results[1];

  outfile.close();
}

template <typename NT> void call_test_benchmark_polytope() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using Hpolytope = HPolytope<Point>;
  {
    Hpolytope P = generate_skinny_cube<Hpolytope>(100, false);
    std::string name = "100_skinny_cube";
    bool centered = false;
    double target_time=20; //secs
    test_benchmark_polytope<NT, Point, Hpolytope>(P, name, false, target_time);
  }

  {
    Hpolytope P = generate_cross<Hpolytope>(5, false);
    std::string name = "5_cross";
    bool centered = false;
    double target_time=10; //secs
    test_benchmark_polytope<NT, Point, Hpolytope>(P, name, centered, target_time);
  }

  {
    Hpolytope P = generate_simplex<Hpolytope>(100, false);
    std::string name = "100_simplex";
    bool centered = false;
    double target_time=15; //secs
    test_benchmark_polytope<NT, Point, Hpolytope>(P, name, centered, target_time);
  }

  {
    Hpolytope P = generate_prod_simplex<Hpolytope>(50, false);
    std::string name = "50_prod_simplex";
    bool centered = false;
    double target_time=15; //secs
    test_benchmark_polytope<NT, Point, Hpolytope>(P, name, centered, target_time);
  }

  {
    Hpolytope P = generate_birkhoff<Hpolytope>(10);
    std::string name = "10_birkhoff";
    bool centered = false;
    double target_time=15; //secs
    test_benchmark_polytope<NT, Point, Hpolytope>(P, name, centered, target_time);
  }

  if (exists_check("netlib/afiro.ine")) {
    Hpolytope P = read_polytope<Hpolytope, NT>("netlib/afiro.ine");
    std::string name = "afiro";
    bool centered = true;
    double target_time=100; //secs
    test_benchmark_polytope<NT, Point, Hpolytope>(P, name, centered, target_time);
  }

  if (exists_check("metabolic_full_dim/polytope_e_coli.ine")) {
    Hpolytope P =
        read_polytope<Hpolytope, NT>("metabolic_full_dim/polytope_e_coli.ine");
    std::string name = "e_coli";
    bool centered = true;
    double target_time=600; //secs
    test_benchmark_polytope<NT, Point, Hpolytope>(P, name, centered, target_time);
  }
}

int main() {

  std::cout
      << "---------------CRHMC polytope sampling benchmarking---------------"
      << std::endl
      << std::endl;
  call_test_benchmark_polytope<double>();
  return 0;
}
