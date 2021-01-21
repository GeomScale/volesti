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

#include "Eigen/Eigen"

#include "ode_solvers/ode_solvers.hpp"
#include "diagnostics/diagnostics.hpp"

#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/known_polytope_generators.h"

template <typename NT>
void run_main() {
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef HPolytope<Point> Hpolytope;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    typedef GaussianFunctor::GradientFunctor<Point> NegativeGradientFunctor;
    typedef GaussianFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
    typedef LeapfrogODESolver<Point, NT, Hpolytope, NegativeGradientFunctor> Solver;
    typedef typename HPolytope<Point>::MT MT;
    typedef typename HPolytope<Point>::VT VT;

    RandomNumberGenerator rng(1);
    unsigned int dim = 10;

    Hpolytope P = generate_birkhoff<Hpolytope>(10);
    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();

    Point x0 = inner_ball.first;
    NT r = inner_ball.second;

    GaussianFunctor::parameters<NT, Point> params(x0, 2 / (r * r), NT(-1));
    GaussianRDHRWalk::Walk<Hpolytope, RandomNumberGenerator> walk(P, x0, params.L, rng);
    int n_warmstart_samples = 0;
    unsigned int walk_length = 30;

    for (int i = 0; i < n_warmstart_samples; i++) {
        walk.apply(P, x0, params.L, walk_length, rng);
    }

    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);

    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);

    // In the first argument put in the address of an H-Polytope
    // for truncated sampling and NULL for untruncated
    HamiltonianMonteCarloWalk::Walk
      <Point, Hpolytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
      hmc(&P, x0, F, f, hmc_params);

    int n_samples = 80000;
    int n_burns = 20000;
    unsigned int min_ess;

    MT samples;
    samples.resize(dim, n_samples - n_burns);

    hmc.solver->eta0 = inner_ball.second / 10;

    for (int i = 0; i < n_samples; i++) {
      if (i % 1000 == 0) std::cerr << ".";
      hmc.apply(rng, 100);
      if (i == n_burns) hmc.disable_adaptive();
      if (i >= n_burns) {
          samples.col(i - n_burns) = hmc.x.getCoefficients();
          std::cout << hmc.x.getCoefficients().transpose() << std::endl;
      }
    }
    std::cerr << std::endl;

    print_diagnostics<NT, VT, MT>(samples, min_ess, std::cerr);

    std::cerr << "Average number of reflections: " <<
        (1.0 * hmc.solver->num_reflections) / hmc.solver->num_steps << std::endl;
    std::cerr << "Step size (final): " << hmc.solver->eta << std::endl;
    std::cerr << "Discard Ratio: " << hmc.discard_ratio << std::endl;
    std::cerr << "Average Acceptance Probability: " << exp(hmc.average_acceptance_log_prob) << std::endl;

}

int main() {
  run_main<double>();
  return 0;
}
