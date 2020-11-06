// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include <string>
#include <typeinfo>
#include <fstream>

#include "doctest.h"
#include "Eigen/Eigen"

#include "ode_solvers.hpp"
#include "diagnostics/geweke.hpp"
#include "diagnostics/multivariate_psrf.hpp"
#include "diagnostics/raftery.hpp"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/known_polytope_generators.h"
#include "misc/misc.h"

struct InnerBallFunctor {

  // Custom density with neg log prob equal to c^T x
  template <
      typename NT,
      typename Point
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number
    Point x0;

    parameters(Point x0_) : order(2), L(1), m(1), kappa(1), x0(x0_) {};

  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT, Point> &params;

    GradientFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * (xs[0] - params.x0);
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      Point y = x - params.x0;
      return 0.5 * y.dot(y);
    }

  };

};

struct CustomFunctor {

  // Custom density with neg log prob equal to || x ||^2 + 1^T x
  template <
      typename NT
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number

    parameters() : order(2), L(2), m(2), kappa(1) {};

    parameters(unsigned int order_) :
      order(order),
      L(2),
      m(2),
      kappa(1)
    {}
  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT> params;

    GradientFunctor() {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-1.0) * Point::all_ones(xs[0].dimension());
        y = y + (-2.0) * xs[0];
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT> params;

    FunctionFunctor() {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      return x.dot(x) + x.sum();
    }

  };

};

template <typename Sampler, typename RandomNumberGenerator, typename NT, typename Point>
void check_ergodic_mean_norm(
    Sampler &sampler,
    RandomNumberGenerator &rng,
    Point &mean,
    unsigned int dim,
    int n_samples=1500,
    int skip_samples=750,
    NT target=NT(0),
    NT tol=5e-1) {

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < n_samples; i++) {
    sampler.apply(rng, 1);
    if (i >= skip_samples) {
      mean = mean + sampler.x;
    }

    #ifdef VOLESTI_DEBUG
      std::cout << sampler.x.getCoefficients().transpose() << std::endl;
    #endif
  }

  auto stop = std::chrono::high_resolution_clock::now();

  long ETA = (long) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  mean = (1.0 / (n_samples - skip_samples)) * mean;

  NT error = abs(NT(mean.dot(mean)) - target);

  if (target != NT(0)) error /= abs(target);

  std::cout << "Dimensionality: " << dim << std::endl;
  std::cout << "Target ergodic mean norm: " << target << std::endl;
  std::cout << "Error (relative if possible) after " << n_samples << " samples: " << error << std::endl;
  std::cout << "ETA (us): " << ETA << std::endl << std::endl;

  CHECK(error < tol);

}

template <typename NT>
void benchmark_hmc(bool truncated) {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;

    NegativeGradientFunctor F;
    NegativeLogprobFunctor f;
    RandomNumberGenerator rng(1);
    unsigned int dim_min = 1;
    unsigned int dim_max = 100;
    int n_samples = 1000;
    long ETA;
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

    for (unsigned int dim = dim_min; dim <= dim_max; dim++) {
      Point x0(dim);
      HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);
      if (truncated) {
        Hpolytope P = generate_cube<Hpolytope>(dim, false);
        HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
        hmc(&P, x0, F, f, hmc_params);
        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n_samples; i++) hmc.apply(rng, 1);
        stop = std::chrono::high_resolution_clock::now();

      } else {

        HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
        hmc(NULL, x0, F, f, hmc_params);

        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n_samples; i++) hmc.apply(rng, 1);
        stop = std::chrono::high_resolution_clock::now();
      }

      ETA = (long) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
      std::cout << ETA << std::endl;
    }

}

template <typename NT>
void test_hmc(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef IsotropicQuadraticFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef IsotropicQuadraticFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;

    IsotropicQuadraticFunctor::parameters<NT> params;
    params.order = 2;

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    RandomNumberGenerator rng(1);
    unsigned int dim = 10;
    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);
    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    Point x0(dim);

    HamiltonianMonteCarloWalk::Walk
      <Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    Point mean(dim);
    check_ergodic_mean_norm(hmc, rng, mean, dim, 75000, 37500, NT(0));
}


template <typename NT>
void test_uld(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef IsotropicQuadraticFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef IsotropicQuadraticFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef RandomizedMipointSDESolver<Point, NT, Hpolytope, NegativeGradientFunctor, RandomNumberGenerator> Solver;

    IsotropicQuadraticFunctor::parameters<NT> params;
    params.order = 2;
    params.alpha = NT(1);

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    RandomNumberGenerator rng(1);
    unsigned int dim = 5;
    UnderdampedLangevinWalk::parameters<NT, NegativeGradientFunctor> uld_params(F, dim);
    Hpolytope P = generate_cube<Hpolytope>(dim, false);
    Point x0(dim);

    UnderdampedLangevinWalk::Walk
      <Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      uld(&P, x0, F, f, uld_params);

    Point mean(dim);
    check_ergodic_mean_norm(uld, rng, mean, dim, 75000, 37500, NT(0));

}

template <typename NT>
void benchmark_netlib() {

    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef boost::mt19937 RNGType;
    typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
    typedef InnerBallFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef InnerBallFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;
    typedef typename HPolytope<Point>::MT MT;
    typedef typename HPolytope<Point>::VT VT;

    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;


    inp.open("./netlib/afiro.ine",std::ifstream::in);
    read_pointset(inp, Pin);
    Hpolytope P(Pin);

    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();

    // Chebyshev center
    Point x0 = inner_ball.first;
    unsigned int dim = x0.dimension();

    InnerBallFunctor::parameters<NT, Point> params(x0);

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    RandomNumberGenerator rng(1);

    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);

    HamiltonianMonteCarloWalk::Walk
      <Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    int n_samples = 80000;
    int n_burns = 0;

    MT samples;
    samples.resize(dim, n_samples - n_burns);

    // hmc.solver->eta0 = 0.5;

    for (int i = 0; i < n_samples; i++) {
      if (i % 1000 == 0) std::cerr << ".";
      hmc.apply(rng, 3);
      if (i >= n_burns) {
          samples.col(i - n_burns) = hmc.x.getCoefficients();
      }
    }
    std::cerr << std::endl;

    std::cout << "Step size (final): " << hmc.solver->eta << std::endl;
    std::cout << "Discard Ratio: " << hmc.discard_ratio << std::endl;
    std::cout << "Average Acceptance Probability: " << exp(hmc.average_acceptance_log_prob) << std::endl;
    std::cout << "PSRF: " <<  multivariate_psrf<NT, VT, MT>(samples) << std::endl;

}

template <typename NT>
void call_test_hmc() {
  std::cout << "--- Testing Hamiltonian Monte Carlo" << std::endl;
  test_hmc<NT>();
}

template <typename NT>
void call_test_uld() {
  std::cout << "--- Testing Underdamped Langevin Diffusion" << std::endl;
  test_uld<NT>();
}

template <typename NT>
void call_test_benchmark_hmc(bool truncated) {
  benchmark_hmc<NT>(truncated);
}

template <typename NT>
void call_test_benchmark_netlib() {
    std::cout << " --- Benchmarking netlib polytopes " << std::endl;
    benchmark_netlib<NT>();
}

TEST_CASE("hmc") {
  call_test_hmc<double>();
}

TEST_CASE("uld") {
  call_test_uld<double>();
}

TEST_CASE("benchmark_hmc") {
  call_test_benchmark_hmc<double>(false);
}

TEST_CASE("benchmark_hmc_truncated") {
  call_test_benchmark_hmc<double>(true);
}

TEST_CASE("benchmark_netlib") {
    call_test_benchmark_netlib<double>();
}
