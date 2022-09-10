// Use backward-mode automatic differentiation using StanMath Library

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou
// Copyright (c) 2022-2022 Zhang zhuyan

// Contributed and/or modified by Zhang zhuyan, as part of Google Summer of Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file


#include <cmath>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include "Eigen/Eigen"

#include "ode_solvers/ode_solvers.hpp"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/known_polytope_generators.h"
#include "diagnostics/diagnostics.hpp"

struct autoDiffFunctor
{

  template <typename T>
  T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta) const
  {

    auto temp = theta.array(); // convert to array so that elementwise operation is possible
    return temp.pow(4).sum() - 2 * temp.pow(2).sum() + 2 * 0.795;
  }
};

struct CustomFunctor
{

  template <
      typename NT>
  struct parameters
  {
    unsigned int order;
    NT L;     // Lipschitz constant for gradient
    NT m;     // Strong convexity constant
    NT kappa; // Condition number
    parameters() : order(2), L(4), m(4), kappa(1){};
  };

  template <
      typename Point>
  struct GradientFunctor
  {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;
    autoDiffFunctor f;
    parameters<NT> &params;

    GradientFunctor(parameters<NT> &params_) : params(params_){};

    // The index i represents the state vector index
    Point operator()(unsigned int const &i, pts const &xs, NT const &t) const
    {
      if (i == params.order - 1)
      {

        double fx;
        Eigen::Matrix<double, Eigen::Dynamic, 1> grad_fx;
        stan::math::gradient(f, xs[0].getCoefficients(), fx, grad_fx);
        Point y(grad_fx * -1);
        return y;
      }
      else
      {
        return xs[i + 1]; // returns derivative
      }
    }
  };

  template <
      typename Point>
  struct FunctionFunctor
  {
    typedef typename Point::FT NT;
    autoDiffFunctor f;

    parameters<NT> &params;

    FunctionFunctor(parameters<NT> &params_) : params(params_){};

    // The index i represents the state vector index
    NT operator()(Point const &x) const
    {
      return f(x.getCoefficients());
    }
  };
};

template <typename NT>
void run_main()
{
  typedef Cartesian<NT> Kernel;
  typedef typename Kernel::Point Point;
  typedef std::vector<Point> pts;
  typedef HPolytope<Point> Hpolytope;
  typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
  typedef CustomFunctor::GradientFunctor<Point> NegativeGradientFunctor;
  typedef CustomFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
  typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;
  typedef typename Hpolytope::MT MT;
  typedef typename Hpolytope::VT VT;
  CustomFunctor::parameters<NT> params;

  NegativeGradientFunctor F(params);
  NegativeLogprobFunctor f(params);

  RandomNumberGenerator rng(1);
  unsigned int dim = 2;

  HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);

  Hpolytope P = generate_cube<Hpolytope>(dim, false);

  Point x0 = -0.25 * Point::all_ones(dim);

  // In the first argument put in the address of an H-Polytope
  // for truncated sampling and NULL for untruncated
  HamiltonianMonteCarloWalk::Walk<Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);
  // hmc(NULL, x0, F, f, hmc_params);

  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
  start = std::chrono::high_resolution_clock::now();

  int n_samples = 50000; // Half will be burned
  int max_actual_draws = n_samples / 2;
  unsigned int min_ess = 0;
  MT samples;
  samples.resize(dim, max_actual_draws);

  for (int i = 0; i < n_samples - max_actual_draws; i++)
  {
    hmc.apply(rng, 3);
  }
  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < max_actual_draws; i++)
  {
    std::cout << hmc.x.getCoefficients().transpose() << std::endl;
    hmc.apply(rng, 3);
    samples.col(i) = hmc.x.getCoefficients();
  }
  stop = std::chrono::high_resolution_clock::now();
  long ETA = (long)std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  std::cerr << "total time taken " << ETA << std::endl;
  std::cerr << std::endl;
  print_diagnostics<NT, VT, MT>(samples, min_ess, std::cerr);
  std::cerr << "min ess " << min_ess << "us" << std::endl;
  std::cerr << "Average time per sample: " << ETA / max_actual_draws << "us" << std::endl;
  std::cerr << "Average time per independent sample: " << ETA / min_ess << "us" << std::endl;
  std::cerr << "Average number of reflections: " << (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps << std::endl;
  std::cerr << "Step size (final): " << hmc.solver->eta << std::endl;
  std::cerr << "Discard Ratio: " << hmc.discard_ratio << std::endl;
  std::cerr << "Average Acceptance Probability: " << exp(hmc.average_acceptance_log_prob) << std::endl;
  std::cerr << "PSRF: " << multivariate_psrf<NT, VT, MT>(samples) << std::endl;
  std::cerr << std::endl;
}

int main()
{
  run_main<double>();
  return 0;
}
