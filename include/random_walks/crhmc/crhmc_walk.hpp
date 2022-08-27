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
#ifndef CRHMC_WALK_HPP
#define CRHMC_WALK_HPP
#include "generators/boost_random_number_generator.hpp"
#include "ode_solvers/ode_solvers.hpp"
#include "random_walks/crhmc/additional_units/auto_tuner.hpp"
#include "random_walks/gaussian_helpers.hpp"
#include <chrono>
struct CRHMCWalk {
  template <typename NT, typename OracleFunctor> struct parameters {
    using Opts = opts<NT>;
    NT epsilon;   // tolerance in mixing
    NT eta = 0.2; // step size
    NT momentum;
    NT effectiveStepSize = 1;
    Opts &options;
    parameters(OracleFunctor const &F, unsigned int dim, Opts &user_options,
               NT epsilon_ = 2)
        : options(user_options) {
      epsilon = epsilon_;
      eta = 1.0 / (dim * sqrt(F.params.L));
      momentum = 1 - std::min(1.0, eta / effectiveStepSize);
    }
  };

  template <typename Point, typename Polytope, typename RandomNumberGenerator,
            typename NegativeGradientFunctor, typename NegativeLogprobFunctor,
            typename Solver>
  struct Walk {
    using point = Point;
    using pts = std::vector<Point>;
    using NT = typename Point::FT;
    using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
    using Sampler = CRHMCWalk::Walk<Point, Polytope, RandomNumberGenerator,
                                    NegativeGradientFunctor,
                                    NegativeLogprobFunctor, Solver>;

    using Opts = typename Polytope::Opts;

    // Hyperparameters of the sampler
    parameters<NT, NegativeGradientFunctor> &params;

    // Numerical ODE solver
    Solver *solver;

    // Dimension
    unsigned int dim;

    // Polytope
    Polytope &P;

    // Discarded Samples
    long total_discarded_samples = 0;
    long num_runs = 0;
    float discard_ratio = 0;

    // Average acceptance probability
    float total_acceptance_prob = 0;
    float average_acceptance_prob = 0;

    // Acceptance probability
    NT prob;
    bool accepted;
    NT accept;
    bool update_modules;

    // References to xs
    Point x, v;

    // Proposal points
    Point x_tilde, v_tilde;

    // Gradient function
    NegativeGradientFunctor &F;

    // Auto tuner
    auto_tuner<Sampler, RandomNumberGenerator> *module_update;

    // Helper variables
    NT H, H_tilde, log_prob, u_logprob;
    // Density exponent
    NegativeLogprobFunctor &f;
#ifdef TIME_KEEPING
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> H_duration =
        std::chrono::duration<double>::zero();
#endif
    Walk(Polytope &Problem, Point &p, NegativeGradientFunctor &neg_grad_f,
         NegativeLogprobFunctor &neg_logprob_f,
         parameters<NT, NegativeGradientFunctor> &param)
        : params(param), F(neg_grad_f), f(neg_logprob_f), P(Problem) {

      dim = p.dimension();

      // Starting point is provided from outside
      x = p;
      accepted = false;
      // Initialize solver
      solver =
          new Solver(0.0, params.eta, pts{x, x}, F, Problem, params.options);
      v = solver->get_state(1);
      module_update = new auto_tuner<Sampler, RandomNumberGenerator>(*this);
      update_modules = params.options.DynamicWeight ||
                       params.options.DynamicRegularizer ||
                       params.options.DynamicStepSize;
    };
    Point GetDirectionWithMomentum(unsigned int const &dim,
                                   RandomNumberGenerator &rng, Point x, Point v,
                                   NT momentum = 0, bool normalize = true) {
      Point z = GetDirection<Point>::apply(dim, rng, normalize);
      solver->ham.move({x, v});
      VT sqrthess = (solver->ham.hess).cwiseSqrt();
      z = Point(sqrthess.cwiseProduct(z.getCoefficients()));
      return v * std::sqrt(momentum) + z * std::sqrt(1 - momentum);
    }
    // Returns the current point in the tranformed in the original space
    inline Point getPoint() { return Point(P.T * x.getCoefficients() + P.y); }
    inline void blendv(Point &x, Point &x_new, std::vector<bool> accepted) {
      for (int i = 0; i < dim; i++) {
        if (accepted[i])
          x(i) = x_new(i);
      }
    }
    inline void apply(RandomNumberGenerator &rng, int walk_length = 1,
                      bool metropolis_filter = true) {

      num_runs++;
      //  Pick a random velocity with momentum
      v = GetDirectionWithMomentum(dim, rng, x, v, params.momentum, false);

      solver->set_state(0, x);
      solver->set_state(1, v);
      // Get proposals
      solver->steps(walk_length, accepted);
      x_tilde = solver->get_state(0);
      v_tilde = solver->get_state(1);

      if (metropolis_filter) {
#ifdef TIME_KEEPING
        start = std::chrono::system_clock::now();
#endif
        // Calculate initial Hamiltonian
        H = solver->ham.hamiltonian(x, v);
        // Calculate new Hamiltonian
        H_tilde = solver->ham.hamiltonian(x_tilde, Point(dim) - v_tilde);
#ifdef TIME_KEEPING
        end = std::chrono::system_clock::now();
        H_duration += end - start;
#endif
        NT feasible = solver->ham.feasible(x_tilde.getCoefficients(),
                                           v_tilde.getCoefficients());
        prob = std::min(1.0, exp(H - H_tilde)) * feasible;

        log_prob = log(prob);
        total_acceptance_prob += prob;

        // Decide to switch
        if (rng.sample_urdist() < prob) {
          x = x_tilde;
          v = v_tilde;
          accepted = true;
        } else {
          total_discarded_samples++;
          accepted = false;
          v = Point(dim) - v;
        }
        discard_ratio = (1.0 * total_discarded_samples) / num_runs;
        average_acceptance_prob = total_acceptance_prob / num_runs;
        accept = accepted ? 1 : 0;
      } else {
        x = x_tilde;
        v = v_tilde;
      }
      if (update_modules) {
        module_update->updateModules(*this, rng);
      }
    }
#ifdef TIME_KEEPING
    void print_timing_information() {
      std::cerr << "--------------Timing Information--------------\n";
      double DU_time = solver->DU_duration.count();
      double DK_time = solver->approxDK_duration.count();
      double H_time = H_duration.count();
      double total_time = H_time + DK_time + DU_time;
      std::cerr << "Total elapsed time: " << total_time << "\n";
      std::cerr << "Computing the Hamiltonian in time, " << H_time << " secs\n";
      std::cerr << "Computing DU partial derivatives in time, " << DU_time
                << " secs\n";
      std::cerr << "Computing DK partial derivatives in time, " << DK_time
                << " secs\n";
    }
#endif
  };
};

#endif
