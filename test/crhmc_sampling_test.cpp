// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer
// of Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file
#include "Eigen/Eigen"
#include "diagnostics/diagnostics.hpp"
#include "doctest.h"
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
#include <atomic>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <thread>
#include <tuple>
#include <typeinfo>
#include <unistd.h>
#include <vector>
#include "preprocess/svd_rounding.hpp"
struct InnerBallFunctor {

  // Gaussian density centered at the inner ball center
  template <typename NT, typename Point>
  struct parameters {
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

  template <typename Point>
  struct GradientFunctor {
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

  template <typename Point>
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_){};

    // The index i represents the state vector index
    NT operator()(Point const &x) const {
      Point y = x - params.x0;
      return 1.0 / (2 * pow(params.sigma, 2)) * y.dot(y);
    }
  };
  template <typename Point>
  struct HessianFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;
    HessianFunctor(parameters<NT, Point> &params_) : params(params_){};

    Point operator()(Point const &x) const {
      return (1.0 / pow(params.sigma, 2)) * Point::all_ones(x.dimension());
    }
  };
};
struct CustomFunctor {

  // Custom density with neg log prob equal to || x ||^2 + 1^T x
  template <typename NT>
  struct parameters {
    unsigned int order;
    NT L;     // Lipschitz constant for gradient
    NT m;     // Strong convexity constant
    NT kappa; // Condition number

    parameters() : order(2), L(2), m(2), kappa(1){};

    parameters(unsigned int order_) : order(order), L(2), m(2), kappa(1) {}
  };

  template <typename Point>
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> params;

    GradientFunctor(){};

    // The index i represents the state vector index
    Point operator()(unsigned int const &i, pts const &xs, NT const &t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * Point::all_ones(xs[0].dimension());
        y = y + (-2.0) * xs[0];
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }
    Point operator()(Point const &x) const {
      Point y = (-1.0) * Point::all_ones(x.dimension());
      y = y + (-2.0) * x;
      return y;
    }
  };

  template <typename Point>
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT> params;

    FunctionFunctor(){};

    // The index i represents the state vector index
    NT operator()(Point const &x) const { return x.dot(x) + x.sum(); }
  };
  template <typename Point>
  struct HessianFunctor {
    typedef typename Point::FT NT;
    Point operator()(Point const &x) const {
      return 2 * Point::all_ones(x.dimension());
    }
  };
};
template <typename NT, typename VT, typename MT>
NT check_interval_psrf(MT &samples, NT target = NT(1.2)) {
  NT max_psrf = NT(0);
  VT intv_psrf = interval_psrf<VT, NT, MT>(samples);
  unsigned int d = intv_psrf.rows();
  for (unsigned int i = 0; i < d; i++) {
    CHECK(intv_psrf(i) < target);
    if (intv_psrf(i) > max_psrf)
      max_psrf = intv_psrf(i);
  }
  return max_psrf;
}

template <typename Sampler, typename RandomNumberGenerator, typename NT,
          typename Point>
void check_ergodic_mean_norm(Sampler &sampler, RandomNumberGenerator &rng,
                             Point &mean, unsigned int dim,
                             int n_samples = 1500, int skip_samples = 750,
                             NT target = NT(0), NT tol = 5e-1) {
  int simdLen = sampler.simdLen;
  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < std::ceil(n_samples / simdLen); i++) {
    sampler.apply(rng, 1);
    if (i >= skip_samples) {
      Point x = Point(sampler.getPoints().rowwise().sum());
      mean = mean + x;
    }

#ifdef VOLESTI_DEBUG
    std::cout << sampler.getPoint().x.getCoefficients().transpose()
              << std::endl;
#endif
  }

  auto stop = std::chrono::high_resolution_clock::now();

  long ETA =
      (long)std::chrono::duration_cast<std::chrono::microseconds>(stop - start)
          .count();

  mean = (1.0 / (n_samples - skip_samples)) * mean;

  NT error = abs(NT(mean.dot(mean)) - target);

  if (target != NT(0))
    error /= abs(target);

  std::cout << "Dimensionality: " << dim << std::endl;
  std::cout << "Target ergodic mean norm: " << target << std::endl;
  std::cout << "Error (relative if possible) after " << n_samples
            << " samples: " << error << std::endl;
  std::cout << "ETA (us): " << ETA << std::endl
            << std::endl;

  CHECK(error < tol);
}
template <typename Polytope, typename NT>
Polytope read_polytope(std::string filename) {
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  inp.open(filename, std::ifstream::in);
  read_pointset(inp, Pin);
  Polytope P(Pin);
  return P;
}

template <typename NT, typename Polytope, int simdLen = 1>
void crhmc_polytope_sampling(
    Polytope &P, NT eta = NT(-1), unsigned int walk_length = 1,
    bool rounding = false, bool centered = false,
    unsigned int max_draws = 80000, unsigned int num_burns = 20000) {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;
  using NegativeGradientFunctor = InnerBallFunctor::GradientFunctor<Point>;
  using NegativeLogprobFunctor = InnerBallFunctor::FunctionFunctor<Point>;
  using HessianFunctor = InnerBallFunctor::HessianFunctor<Point>;
  using MT = typename Polytope::MT;
  using VT = typename Polytope::VT;
  using Input = crhmc_input<MT, Point, NegativeLogprobFunctor,
                            NegativeGradientFunctor, HessianFunctor>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using Opts = opts<NT>;
  using Solver = ImplicitMidpointODESolver<Point, NT, CrhmcProblem,
                                           NegativeGradientFunctor, simdLen>;

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

  NegativeGradientFunctor F(params);
  NegativeLogprobFunctor f(params);
  HessianFunctor H(params);
  int max_actual_draws = max_draws - num_burns;
  unsigned int min_ess = 0;

  MT samples;
  samples.resize(dim, max_actual_draws);
  NT ETA;
  NT max_psrf;

  Opts options;
  options.simdLen = simdLen;
  CRHMCWalk::parameters<NT, NegativeGradientFunctor> crhmc_params(F, dim,
                                                                  options);
  Input input = Input(P.dimension(), f, F, H);
  input.Aineq = P.get_mat();
  input.bineq = P.get_vec();

  CrhmcProblem crhmc_problem = CrhmcProblem(input);
  Point x_start(crhmc_problem.center);
  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator,
                  NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
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

  for (unsigned int i = 0; i < std::ceil(max_actual_draws / simdLen); i++) {
    for (int k = 0; k < walk_length; k++) {
      crhmc.apply(rng, 1);
    }
    MT sample = crhmc.getPoints();
    if (i * simdLen + simdLen - 1 < max_actual_draws) {
      samples(Eigen::all, Eigen::seq(i * simdLen, i * simdLen + simdLen - 1)) = sample;
    } else {
      samples(Eigen::all, Eigen::seq(i * simdLen, max_actual_draws - 1)) = sample(Eigen::all, Eigen::seq(0, max_actual_draws - 1 - simdLen * i));
    }
    if (i % 1000 == 0 && i > 0)
      std::cout << ".";
  }
  std::cout << std::endl;
  std::cout << "Step size (final): " << crhmc.solver->eta << std::endl;
  std::cout << "Discard Ratio: " << crhmc.discard_ratio << std::endl;
  std::cout << "Average Acceptance Probability: "
            << crhmc.average_acceptance_prob << std::endl;
  max_psrf = check_interval_psrf<NT, VT, MT>(samples);
  std::cout << "max_psrf: " << max_psrf << std::endl;
  std::cout << std::endl;
}
inline bool exists_check(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

template <typename NT, typename Point, typename HPolytope, int simdLen = 1>
void test_sampling_polytope(HPolytope &P, std::string &name, bool centered,
                            int walk_length = 1) {
  NT step_size = 0;
  std::pair<Point, NT> inner_ball;
  std::cout << name << std::endl;
  P.normalize();
  inner_ball = P.ComputeInnerBall();
  step_size = inner_ball.second / 10;
  crhmc_polytope_sampling<NT, HPolytope, simdLen>(P, step_size, walk_length, false, centered);
}
template <typename NT, int simdLen = 1>
void call_test_sampling_polytope() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using Hpolytope = HPolytope<Point>;
  std::cout << " ---Sampling polytopes " << std::endl;

  {
    Hpolytope P = generate_skinny_cube<Hpolytope>(100, false);
    std::string name = "100_skinny_cube";
    bool centered = false;
    test_sampling_polytope<NT, Point, Hpolytope, simdLen>(P, name, false);
  }

  {
    Hpolytope P = generate_cross<Hpolytope>(5, false);
    std::string name = "5_cross";
    bool centered = false;
    test_sampling_polytope<NT, Point, Hpolytope, simdLen>(P, name, centered);
  }

  {
    Hpolytope P = generate_simplex<Hpolytope>(100, false);
    std::string name = "100_simplex";
    bool centered = false;
    test_sampling_polytope<NT, Point, Hpolytope, simdLen>(P, name, centered);
  }

  {
    Hpolytope P = generate_prod_simplex<Hpolytope>(50, false);
    std::string name = "50_prod_simplex";
    bool centered = false;
    test_sampling_polytope<NT, Point, Hpolytope, simdLen>(P, name, centered);
  }

  {
    Hpolytope P = generate_birkhoff<Hpolytope>(10);
    std::string name = "10_birkhoff";
    bool centered = false;
    test_sampling_polytope<NT, Point, Hpolytope, simdLen>(P, name, centered);
  }
}

template <typename NT>
void benchmark_cube_crhmc() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using pts = std::vector<Point>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using NegativeGradientFunctor = CustomFunctor::GradientFunctor<Point>;
  using NegativeLogprobFunctor = CustomFunctor::FunctionFunctor<Point>;
  using Input =
      crhmc_input<MT, Point, NegativeLogprobFunctor, NegativeGradientFunctor>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;
  using Solver = ImplicitMidpointODESolver<Point, NT, CrhmcProblem,
                                           NegativeGradientFunctor>;
  using Opts = opts<NT>;
  NegativeGradientFunctor g;
  NegativeLogprobFunctor f;
  RandomNumberGenerator rng(1);
  Opts options;
  unsigned int dim_min = 1;
  unsigned int dim_max = 100;
  int n_samples = 1000;
  long ETA;
  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

  for (unsigned int dim = dim_min; dim <= dim_max; dim++) {
    CRHMCWalk::parameters<NT, NegativeGradientFunctor> crhmc_params(g, dim,
                                                                    options);
    Input input = Input(dim, f, g);
    input.lb = -VT::Ones(dim);
    input.ub = VT::Ones(dim);
    CrhmcProblem P = CrhmcProblem(input);
    Point x0(P.center);
    CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator,
                    NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
        crhmc(P, x0, g, f, crhmc_params);
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_samples; i++)
      crhmc.apply(rng, 1);
    stop = std::chrono::high_resolution_clock::now();

    ETA = (long)std::chrono::duration_cast<std::chrono::microseconds>(stop -
                                                                      start)
              .count();
    std::cout << ETA << std::endl;
  }
}

template <typename NT, int simdLen = 1>
void test_crhmc() {
  using Kernel = Cartesian<NT>;
  using Point = typename Kernel::Point;
  using pts = std::vector<Point>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using NegativeGradientFunctor =
      IsotropicQuadraticFunctor::GradientFunctor<Point>;
  using NegativeLogprobFunctor =
      IsotropicQuadraticFunctor::FunctionFunctor<Point>;
  using Input =
      crhmc_input<MT, Point, NegativeLogprobFunctor, NegativeGradientFunctor>;
  using CrhmcProblem = crhmc_problem<Point, Input>;
  using RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, NT>;
  using Solver = ImplicitMidpointODESolver<Point, NT, CrhmcProblem, NegativeGradientFunctor, simdLen>;
  using Opts = opts<NT>;
  IsotropicQuadraticFunctor::parameters<NT> params;
  params.order = 2;
  NegativeGradientFunctor g(params);
  NegativeLogprobFunctor f(params);
  RandomNumberGenerator rng(1);
  unsigned int dim = 10;
  Opts options;
  options.simdLen = simdLen;

  CRHMCWalk::parameters<NT, NegativeGradientFunctor> crhmc_params(g, dim,
                                                                  options);
  Input input = Input(dim, f, g);
  input.lb = -VT::Ones(dim);
  input.ub = VT::Ones(dim);
  CrhmcProblem P = CrhmcProblem(input);
  Point x0(dim);
  CRHMCWalk::Walk<Point, CrhmcProblem, RandomNumberGenerator,
                  NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      crhmc(P, x0, g, f, crhmc_params);
  Point mean(dim);
  check_ergodic_mean_norm(crhmc, rng, mean, dim, 75000, 37500, NT(0));
}

template <typename NT>
void call_test_crhmc() {
  std::cout << "--- Testing Constrained Riemannian Hamiltonian Monte Carlo"
            << std::endl;
  std::cout << "------------SIMDLEN=1-------------------\n"
            << std::endl;
  test_crhmc<NT, 1>();
  std::cout << "------------SIMDLEN=4-------------------\n"
            << std::endl;
  test_crhmc<NT, 4>();
}
template <typename NT>
void call_test_benchmark_cube_crhmc() {
  benchmark_cube_crhmc<NT>();
}

TEST_CASE("crhmc") {
  call_test_crhmc<double>();
}

TEST_CASE("benchmark_crhmc_cube") {
  call_test_benchmark_cube_crhmc<double>();
}

TEST_CASE("test_polytope_sampling_crhmc") {
  std::cout << "------------SIMDLEN=1-------------------\n"
            << std::endl;
  call_test_sampling_polytope<double, 1>();
  std::cout << "------------SIMDLEN=4-------------------\n"
            << std::endl;
  call_test_sampling_polytope<double, 4>();
}
