// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Lee, Yin Tat, Ruoqi Shen, and Kevin Tian. "Logsmooth Gradient Concentration
// and Tighter Runtimes for Metropolized Hamiltonian Monte Carlo."
#ifndef HAMILTONIAN_MONTE_CARLO_WALK_HPP
#define HAMILTONIAN_MONTE_CARLO_WALK_HPP


#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"
#include "ode_solvers/ode_solvers.hpp"

struct HamiltonianMonteCarloWalk {

  template
  <
    typename NT,
    typename OracleFunctor
  >
  struct parameters {
    NT epsilon; // tolerance in mixing
    NT eta; // step size

    parameters(
      OracleFunctor const& F,
      unsigned int dim,
      NT epsilon_=2)
    {
      epsilon = epsilon_;
      eta = 1.0 / (dim * sqrt(F.params.L));
      // eta = 1.0 /
      //   (sqrt(20 * F.params.L * pow(dim, 3)));
    }
  };

  template
  <
    typename Point,
    typename Polytope,
    typename RandomNumberGenerator,
    typename NegativeGradientFunctor,
    typename NegativeLogprobFunctor,
    typename Solver
  >
  struct Walk {

    typedef std::vector<Point> pts;
    typedef typename Point::FT NT;
    typedef std::vector<Polytope*> bounds;

    // Hyperparameters of the sampler
    parameters<NT, NegativeGradientFunctor> &params;

    // Numerical ODE solver
    Solver *solver;

    // Dimension
    unsigned int dim;

    // Discarded Samples
    long total_discarded_samples = 0;
    long num_runs = 0;
    float discard_ratio = 0;

    // Average acceptance probability
    float total_acceptance_log_prob = 0;
    float average_acceptance_log_prob = 0;

    // References to xs
    Point x, v;

    // Proposal points
    Point x_tilde, v_tilde;

    // Gradient function
    NegativeGradientFunctor &F;

    bool accepted;

    // Helper variables
    NT H, H_tilde, log_prob, u_logprob;

    // Density exponent
    NegativeLogprobFunctor &f;

    Walk(Polytope *P,
      Point &p,
      NegativeGradientFunctor &neg_grad_f,
      NegativeLogprobFunctor &neg_logprob_f,
      parameters<NT, NegativeGradientFunctor> &param) :
      params(param),
      F(neg_grad_f),
      f(neg_logprob_f)
    {

      dim = p.dimension();

      // Starting point is provided from outside
      x = p;

      accepted = false;

      // Initialize solver
      solver = new Solver(0, params.eta, pts{x, x}, F, bounds{P, NULL});

    };

    inline void apply(RandomNumberGenerator &rng,
                      int walk_length=1,
                      bool metropolis_filter=true)
    {
      num_runs++;

      // Pick a random velocity
      v = GetDirection<Point>::apply(dim, rng, false);

      solver->set_state(0, x);
      solver->set_state(1, v);

      // Get proposals
      solver->steps(walk_length, accepted);
      x_tilde = solver->get_state(0);
      v_tilde = solver->get_state(1);

      if (metropolis_filter) {
        // Calculate initial Hamiltonian
        H = hamiltonian(x, v);

        // Calculate new Hamiltonian
        H_tilde = hamiltonian(x_tilde, v_tilde);

        // Log-sum-exp trick
        log_prob = H - H_tilde < 0 ? H - H_tilde : 0;

        // Decide to switch
        u_logprob = log(rng.sample_urdist());
        total_acceptance_log_prob += log_prob;
        if (u_logprob < log_prob) {
          x = x_tilde;
          accepted = true;
        }
        else {
          total_discarded_samples++;
          accepted = false;
        }
      } else {
        x = x_tilde;
        accepted = true;
      }

      discard_ratio = (1.0 * total_discarded_samples) / num_runs;
      average_acceptance_log_prob = total_acceptance_log_prob / num_runs;

    }

    inline NT hamiltonian(Point &pos, Point &vel) const {
      return f(pos) + 0.5 * vel.dot(vel);
    }

    void disable_adaptive() {
      solver->disable_adaptive();
    }

    void enable_adaptive() {
      solver->enable_adaptive();
    }
  };
};

#endif // HAMILTONIAN_MONTE_CARLO_WALK_HPP
