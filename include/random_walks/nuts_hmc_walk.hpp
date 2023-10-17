// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2022 Vissarion Fisikopoulos
// Copyright (c) 2018-2022 Apostolos Chalkis
// Copyright (c) 2020-2022 Elias Tsigaridas
// Copyright (c) 2020-2022 Marios Papachristou

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Matthew D. Hoffman, Andrew Gelman. "The No-U-Turn Sampler: 
// Adaptively Setting Path Lengths in Hamiltonian Monte Carlo", 2011.

// Comment: Compared to [Matthew D. Hoffman, Andrew Gelman, 2011]
// we modify the step of Nesterov's algorithm in the burn in phase.

#ifndef NUTS_HAMILTONIAN_MONTE_CARLO_WALK_HPP
#define NUTS_HAMILTONIAN_MONTE_CARLO_WALK_HPP


#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"
#include "ode_solvers/ode_solvers.hpp"
#include "preprocess/estimate_L_smooth_parameter.hpp"

struct NutsHamiltonianMonteCarloWalk {

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
      eta = F.params.L > 0 ? 10.0 / (dim * sqrt(F.params.L)) : 0.005;
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
    long num_runs = 0;
    long total_acceptance = 0;

    // Average acceptance probability
    NT average_acceptance = 0;

    // References to xs
    Point x, v;

    // Helper points
    Point v_pl, v_min, v_min_j, v_pl_j, X_pl, X_pl_j, X_min, X, X_rnd_j, X_min_j, x_pl_min;

    // Gradient function
    NegativeGradientFunctor &F;

    bool accepted;

    // Burnin parameters
    NT eps_step, mu, log_tilde_eps, H_tilde, alpha, na;
    const NT delta = NT(0.65), Delta_max = NT(1000), gamma = NT(0.05), t0 = NT(10), kk = NT(0.85);

    // Density exponent
    NegativeLogprobFunctor &f;

    Walk(Polytope *P,
         Point &p,
         NegativeGradientFunctor &neg_grad_f,
         NegativeLogprobFunctor &neg_logprob_f,
         parameters<NT, NegativeGradientFunctor> &param,
         bool burn_in_phase = true) :
         params(param),
         F(neg_grad_f),
         f(neg_logprob_f)
    {
      dim = p.dimension();

      v_pl.set_dimension(dim);
      v_min.set_dimension(dim);
      v_min_j.set_dimension(dim);
      v_pl_j.set_dimension(dim);
      X_pl.set_dimension(dim);
      X_pl_j.set_dimension(dim);
      X_min.set_dimension(dim);
      X.set_dimension(dim);
      X_rnd_j.set_dimension(dim);
      X_min_j.set_dimension(dim);
      x_pl_min.set_dimension(dim);

      eps_step = params.eta;
      mu = std::log(10*eps_step);
      log_tilde_eps = NT(0);
      H_tilde = NT(0);
      alpha = NT(0);
      na = NT(0);

      // Starting point is provided from outside
      x = p;

      accepted = false;

      // Initialize solver
      solver = new Solver(0, params.eta, pts{x, x}, F, bounds{P, NULL});
      disable_adaptive();

      if (burn_in_phase)
      {
        RandomNumberGenerator rng(dim);
        burnin(rng);
      }
    };


    inline void burnin(RandomNumberGenerator &rng,
                       unsigned int N = 1000,
                       unsigned int walk_length=1)
    {
      reset_num_runs();
      Point p = x;
      NT L;

      if ((solver->get_bounds())[0] == NULL) 
      {
        L = (NT(100) / NT(dim)) * (NT(100) / NT(dim));
      }
      else
      {
        Polytope K = *(solver->get_bounds())[0];
        L = estimate_L_smooth(K, p, walk_length, F, rng);
      }

      eps_step = NT(5) / (NT(dim) * std::sqrt(L));
      solver->set_eta(eps_step);

      for (int i = 0; i < N; i++)
      {
        apply(rng, walk_length, true);
        solver->set_eta(eps_step);
      }
      reset_num_runs();
    }


    inline void apply(RandomNumberGenerator &rng,
                      unsigned int walk_length=1,
                      bool burnin = false)
    {
      num_runs++;

      int x_counting_total = 0;

      // Pick a random velocity
      v = GetDirection<Point>::apply(dim, rng, false);
      
      v_pl = v;
      v_min = NT(-1) * v;
      X_pl = x;
      X_min = x;
      
      NT h1 = hamiltonian(x, v);

      NT uu = std::log(rng.sample_urdist()) - h1;
      int j = -1;
      bool s = true;
      bool updated = false;
      bool pos_state_single = false;

      if (burnin)
      {
        alpha = NT(0);
      }

      while (s)
      {
        j++;
        
        if (burnin)
        {
          na = std::pow(NT(2), NT(j));
        }

        NT dir = rng.sample_urdist();
                
        if (dir > 0.5)
        {
          v = v_pl;
          X = X_pl;
        }
        else
        {
          v = v_min;
          X = X_min;
        }
        X_rnd_j = X;
                
        int x_counting = 0;
        int num_samples = int(std::pow(NT(2), NT(j)));
        accepted = false;

        for (int k = 1; k <= num_samples; k++)
        {
          if (!accepted)
          {
            solver->set_state(0, X);
            solver->set_state(1, v);
          }

          // Get proposals
          solver->steps(walk_length, accepted);
          accepted = true;

          X = solver->get_state(0);
          v = solver->get_state(1);

          NT hj = hamiltonian(X, v);

          if (burnin)
          {
            alpha += std::min(NT(1), std::exp(-hj + h1));
          }

          if (uu > Delta_max - hj)
          {
            s = false;
            break;
          }

          bool pos_state = false;
          if (uu < -hj) 
          {
            pos_state = true;
            pos_state_single = true;
            x_counting = x_counting + 1;
            x_counting_total = x_counting_total + 1;
          }
                 
          if (k == 1) 
          {
            if (dir > 0.5) 
            {
              X_min_j = X;
              v_min_j = v;
            } 
            else 
            {
              X_pl_j = X;
              v_pl_j = v;
            }
          }
          if (k == num_samples) 
          {
            if (dir > 0.5) 
            {
              x_pl_min = X - X_min_j;
              if ((x_pl_min.dot(v) < 0) || (x_pl_min.dot(v_min_j) < 0)) 
              {
                s = false;
              }
            } 
            else 
            {
              x_pl_min = X_pl_j - X;
              if ((x_pl_min.dot(v) < 0) || (x_pl_min.dot(v_pl_j) < 0)) 
              {
                s = false;
              }
            }
          }
          if ((rng.sample_urdist() < (1/NT(x_counting))) && pos_state) 
          {
            X_rnd_j = X;
          }
        }

        if (dir > 0.5) 
        {
          X_pl = X;
          v_pl = v;
        } 
        else 
        {
          X_min = X;
          v_min = v;
        }
                
        if (s && (rng.sample_urdist() < (NT(x_counting) / NT(x_counting_total)))) 
        {
          x = X_rnd_j;
          if (pos_state_single)
          { 
            updated = true;
          }
        }
                
        if (s) 
        {
          x_pl_min = X_pl - X_min;
          if ((x_pl_min.dot(v_min) < 0) || (x_pl_min.dot(v_pl) < 0)) 
          {
              s = false;
          }
        }
      }

      if (updated)
      {
        total_acceptance++;
      }
      average_acceptance = NT(total_acceptance) / NT(num_runs);

      if (burnin)
      {
        H_tilde = (NT(1) - NT(1) / (NT(num_runs) + t0)) * H_tilde + (NT(1) / (NT(num_runs) + t0)) * (delta - alpha / na);
        NT log_eps = mu - (std::sqrt(NT(num_runs)) / gamma) * H_tilde;

        // TODO: use the following to generalize Nesterov's algorithm
        //log_tilde_eps = std::pow(mu,-kk) * log_eps + (NT(1) - std::pow(mu,-kk))*log_tilde_eps;

        eps_step = std::exp(log_eps);
      }
    }

    inline NT hamiltonian(Point &pos, Point &vel) const {
      return f(pos) + 0.5 * vel.dot(vel);
    }

    inline NT get_eta_solver() {
      return solver->get_eta();
    }

    void disable_adaptive() {
      solver->disable_adaptive();
    }

    void enable_adaptive() {
      solver->enable_adaptive();
    }

    void reset_num_runs() {
      num_runs = 0;
      total_acceptance = 0;
    }

    NT get_ratio_acceptance() {
      return average_acceptance;
    }
  };
};

#endif // HAMILTONIAN_MONTE_CARLO_WALK_HPP
