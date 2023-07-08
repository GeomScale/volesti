// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Shen, Ruoqi, and Yin Tat Lee. "The randomized midpoint method for
// log-concave sampling." Advances in Neural Information Processing Systems. 2019.

#ifndef LANGEVIN_WALK_HPP
#define LANGEVIN_WALK_HPP

#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"

struct UnderdampedLangevinWalk {

  template
  <
    typename NT,
    typename OracleFunctor
  >
  struct parameters {
    NT epsilon; // tolerance in mixing
    NT eta; // step size
    NT u;
    parameters(
      OracleFunctor const& F,
      unsigned int dim,
      NT epsilon_=1e-4)
    {
      epsilon = epsilon_;
      u = 1.0 / F.params.L;
      eta = 1.0 / (sqrt(20 * F.params.L));

      // eta = std::min(pow(epsilon, 1.0 / 3) /
      //                       pow(F.params.kappa, 1.0 / 6) *
      //                       pow(log(1.0 / epsilon), - 1.0 / 6),
      //                       pow(epsilon, 2.0 / 3) *
      //                       pow(log(1.0 / epsilon), - 1.0 / 3));
    }
  };

  template
  <
    typename Point,
    typename Polytope,
    typename RandomNumberGenerator,
    typename NegativeGradientFunctor,
    typename NegativeLogprobFunctor,
    typename Solver_ // Dummy template argument
  >
  struct Walk {

    typedef std::vector<Point> pts;
    typedef typename Point::FT NT;
    typedef std::vector<Polytope*> bounds;
    typedef RandomizedMipointSDESolver<Point, NT, Polytope, NegativeGradientFunctor, RandomNumberGenerator> Solver;

    parameters<NT, NegativeGradientFunctor> &params;

    // Numerical ODE solver
    Solver *solver;

    unsigned int dim;

    // References to xs
    Point x, v;
    // Proposal points
    Point x_tilde, v_tilde;

    // Gradient Function
    NegativeGradientFunctor &F;

    // Density exponent
    NegativeLogprobFunctor &f;

    NT H_tilde, H, log_prob, u_logprob;

    Walk(Polytope *P,
      Point &initial_x,
      NegativeGradientFunctor &neg_grad_f,
      NegativeLogprobFunctor &neg_logprob_f,
      parameters<NT, NegativeGradientFunctor> &param) :
      params(param),
      F(neg_grad_f),
      f(neg_logprob_f)
    {
      // Starting point is provided from outside
      x = initial_x;
      dim = initial_x.dimension();
      v = Point(dim);

      solver = new Solver(0, params.eta, pts{x, v}, F, bounds{P, NULL}, params.u);

    };

    inline void apply(
      RandomNumberGenerator &rng,
      int walk_length=1,
      bool metropolis_filter=false)
    {
      solver->set_state(0, x);
      solver->set_state(1, v);

      // Get proposals
      solver->steps(walk_length, rng);
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
        if (u_logprob < log_prob) {
          x = x_tilde;
          v = v_tilde;
        }
      } else {
        x = x_tilde;
        v = v_tilde;
      }

    }

    inline NT hamiltonian(Point &pos, Point &vel) const {
      return f(pos) + 1.0 / (2 * params.u) * vel.dot(vel);
    }

    void disable_adaptive() {
        // TODO Implement
    }

    void enable_adaptive() {
        // TODO Implement
    }
  };
};

#endif // LANGEVIN_WALK_HPP
