// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2022 Vissarion Fisikopoulos
// Copyright (c) 2018-2022 Apostolos Chalkis
// Copyright (c) 2020-2022 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file

// References

#ifndef NUTS_HAMILTONIAN_MONTE_CARLO_WALK_HPP
#define NUTS_HAMILTONIAN_MONTE_CARLO_WALK_HPP


#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_helpers.hpp"
#include "ode_solvers/ode_solvers.hpp"

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
    Point v_pl, v_min, v_min_j, v_pl_j, X_pl, X_pl_j, X_min, X, X_rnd_j, X_min_j, x_pl_min;

    // Gradient function
    NegativeGradientFunctor &F;

    bool accepted;

    // Helper variables
    NT H, H_tilde, log_prob, u_logprob, Delta_max;

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

      Delta_max = NT(1000)

      // Starting point is provided from outside
      x = p;

      accepted = false;

      // Initialize solver
      solver = new Solver(0, params.eta, pts{x, x}, F, bounds{P, NULL});

    };

    inline void apply(
      RandomNumberGenerator &rng,
      int walk_length=1,
      bool metropolis_filter=true)
    {
      num_runs++;

      int x_counting_total = 0;

      // Pick a random velocity
      v = GetDirection<Point>::apply(dim, rng, false);
      
      v_pl = v;
      v_min = -v;
      X_pl = x;
      X_min = x;
      
      //grad_x_pl = esti_grad_invW_opt(f_utils, X0, T0i);
      //grad_x_min = grad_x_pl;
      NT h1 = hamiltonian(x,v)
      //h1 = f(X0, T0i, f_utils) + 0.5 * (v * v);

      NT uu = log(rng.sample_urdist()) - h1;
      int j = -1;
      int s = 1;

      while (s == 1)
      {
        j++;// j + 1;
        dir = rng.sample_urdist();
        //%v = v*dir;
                
        if (dir > 0.5)
        {
          //grad_x = grad_x_pl;
          v = v_pl;
          X = X_pl;
          //Ti = Ti_pl;
        }
        else
        {
          //grad_x = grad_x_min;
          v = v_min;
          X = X_min;
          //Ti = Ti_min;
        }
        X_rnd_j = X;
        //Ti_rnd_j = Ti;
        //%v_rnd_j = v;
                
        int x_counting = 0;
        int num_samples = int(std::pow(2, j));
        for (int k = 1; k <= num_samples; k++)// = 1:2^j)
        {
          //v = v - (eta/2) * grad_x;
          //T = eta;

          solver->set_state(0, X);
          solver->set_state(1, v);

          // Get proposals
          solver->steps(walk_length, accepted);
          X = solver->get_state(0);
          v = solver->get_state(1);

          //hj = f(X, Ti, f_utils) + 0.5 * (v * v);
          NT hj = hamiltonian(X,v);
          if (uu > Delta_max - hj){
            s = 0;
            break;
          }
          bool pos_state = false;
          if (uu < -hj) {
            pos_state = true;
            x_counting = x_counting + 1;
            x_counting_total = x_counting_total + 1;
          }
          //%pos_state
                 
          if (k==1) {
            if (dir > 0.5) {
              X_min_j = X;
              //Ti_min_j = Ti;
              v_min_j = v;
            } else {
              X_pl_j = X;
              //Ti_pl_j = Ti;
              v_pl_j = v;
            }
          }
          if (k == num_samples) {
            if (dir > 0.5) {
              //x_pl = [X(lower); Ti];
              //x_min = [X_min_j(lower); Ti_min_j];
              x_pl_min = X - X_min_j;
              if ((x_pl_min.dot(v) < 0) || (x_pl_min.dot(v_min_j) < 0)) {
                s = 0;
              }
            } else {
              //x_pl = [X_pl_j(lower); Ti_pl_j];
              //x_min = [X(lower); Ti];
              x_pl_min = X_pl_j - X;
              if ((x_pl_min.dot(v) < 0) || (x_pl_min.dot(v_pl_j) < 0)) {
                s = 0;
              }
            }
          }
          if (rng.sample_urdist() < (1/NT(x_counting)) && pos_state) {
            X_rnd_j = X;
            //Ti_rnd_j = Ti;
            //%v_rnd_j = v;
          }
        }

        if (dir > 0.5) {
          X_pl = X;
          v_pl = v;
          //Ti_pl = Ti;
          //grad_x_pl = grad_x;
        } else {
          X_min = X;
          v_min = v;
          //Ti_min = Ti;
          //grad_x_min = grad_x;
        }
                
        if (s == 1 && rng.sample_urdist() < (NT(x_counting) / NT(x_counting_total))) {
          x = X_rnd_j;
          //T0i = Ti_rnd_j;
          accepted = true;
          //%v_nx = v_rnd_j;
        }
                
        if (s==1) {
          //x_pl = [X_pl(lower); Ti_pl];
          //x_min = [X_min(lower); Ti_min];
          //%(x_pl - x_min)*v_min
          //%(x_pl - x_min)*v_pl
          x_pl_min = X_pl - X_min;
          if ((x_pl_min.dot(v_min) < 0) || (x_pl_min.dot(v_pl) < 0)) {
              s = 0;
          }
        }
      }

      //x = X0;
      
      /*if (metropolis_filter) {
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
      }*/

      //discard_ratio = (1.0 * total_discarded_samples) / num_runs;
      //average_acceptance_log_prob = total_acceptance_log_prob / num_runs;

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
