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
  template
  <
    typename NT,
    typename OracleFunctor
  >
  struct parameters {
    using Opts = opts<NT>;
    NT epsilon;   // tolerance in mixing
    NT eta = 0.2; // step size
    NT momentum;
    NT effectiveStepSize = 1;
    Opts &options;
    parameters(OracleFunctor const &F,
      unsigned int dim,
      Opts &user_options,
      NT epsilon_ = 2)  :
      options(user_options)
    {
      epsilon = epsilon_;
      eta = 1.0 / (dim * sqrt(F.params.L));
      momentum = 1 - std::min(1.0, eta / effectiveStepSize);
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
    using point = Point;
    using pts = std::vector<Point>;
    using NT = typename Point::FT;
    using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
    using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using Sampler = CRHMCWalk::Walk<Point, Polytope, RandomNumberGenerator,
                                    NegativeGradientFunctor,
                                    NegativeLogprobFunctor, Solver>;

    using Opts = typename Polytope::Opts;
    using IVT = Eigen::Matrix<int, Eigen::Dynamic, 1>;

    // Hyperparameters of the sampler
    parameters<NT, NegativeGradientFunctor> &params;

    // Numerical ODE solver
    std::unique_ptr<Solver> solver;

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
    VT prob;
    bool accepted;
    IVT accept;
    bool update_modules;
    int simdLen;
    // References to xs
    MT x, v;

    // Proposal points
    MT x_tilde, v_tilde;

    // Gradient function
    NegativeGradientFunctor &F;

    // Auto tuner
    std::unique_ptr<auto_tuner<Sampler, RandomNumberGenerator>>module_update;

    // Helper variables
    VT H, H_tilde;
    // Density exponent
    NegativeLogprobFunctor &f;
#ifdef TIME_KEEPING
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> H_duration = std::chrono::duration<double>::zero();
#endif
    Walk(Polytope &Problem,
      Point &p,
      NegativeGradientFunctor &neg_grad_f,
      NegativeLogprobFunctor &neg_logprob_f,
      parameters<NT, NegativeGradientFunctor> &param) :
      params(param),
      P(Problem),
      F(neg_grad_f),
      f(neg_logprob_f)
    {

      dim = p.dimension();
      simdLen = params.options.simdLen;
      // Starting point is provided from outside
      x = p.getCoefficients() * MT::Ones(1, simdLen);
      accepted = false;
      // Initialize solver
      solver = std::unique_ptr<Solver>(new Solver(0.0, params.eta, {x, x}, F, Problem, params.options));
      v = MT::Zero(dim, simdLen);
      module_update = std::unique_ptr<auto_tuner<Sampler, RandomNumberGenerator>>(new auto_tuner<Sampler, RandomNumberGenerator>(*this));
      update_modules = params.options.DynamicWeight ||
                       params.options.DynamicRegularizer ||
                       params.options.DynamicStepSize;
    };
    // Sample a new velocity with momentum
    MT get_direction_with_momentum(unsigned int const &dim,
                                RandomNumberGenerator &rng, MT const &x, MT v,
                                NT momentum = 0, bool normalize = true)
    {
      MT z = MT(dim, simdLen);
      for (int i = 0; i < simdLen; i++)
      {
        z.col(i) = GetDirection<Point>::apply(dim, rng, normalize).getCoefficients();
      }
      solver->ham.move({x, v});
      MT sqrthess = (solver->ham.hess).cwiseSqrt();
      z = sqrthess.cwiseProduct(z);
      return v * std::sqrt(momentum) + z * std::sqrt(1 - momentum);
    }
    // Returns the current point in the tranformed in the original space
    inline MT getPoints() { return (P.T * x).colwise() + P.y; }
    // Returns the current point in the tranformed in the original space
    inline Point getPoint() { return Point(P.T * x.col(0) + P.y); }

    inline MT masked_choose(MT &x, MT &x_tilde, IVT &accept) {
      return accept.transpose().replicate(x.rows(), 1).select(x_tilde, x);
    }
    inline void disable_adaptive(){
      update_modules=false;
    }
    inline void apply(RandomNumberGenerator &rng,
      int walk_length = 1,
      bool metropolis_filter = true)
    {
      num_runs++;
      //  Pick a random velocity with momentum
      v = get_direction_with_momentum(dim, rng, x, v, params.momentum, false);
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
        H_tilde = solver->ham.hamiltonian(x_tilde, -v_tilde);

#ifdef TIME_KEEPING
        end = std::chrono::system_clock::now();
        H_duration += end - start;
#endif
        VT feasible = solver->ham.feasible(x_tilde,
                                           v_tilde);
        prob = (1.0 < exp((H - H_tilde).array())).select(1.0, exp((H - H_tilde).array()));
        prob = (feasible.array() > 0.5).select(prob, 0);

        total_acceptance_prob += prob.sum();
        VT rng_vector = VT(simdLen);
        for (int i = 0; i < simdLen; i++)
        {
          rng_vector(i) = rng.sample_urdist();
        }
        accept = (rng_vector.array() < prob.array()).select(1 * IVT::Ones(simdLen), 0 * IVT::Ones(simdLen));

        x = masked_choose(x, x_tilde, accept);
        v = -v;
        v = masked_choose(v, v_tilde, accept);
        total_discarded_samples += simdLen - accept.sum();
        discard_ratio = (1.0 * total_discarded_samples) / (num_runs * simdLen);
        average_acceptance_prob = total_acceptance_prob / (num_runs * simdLen);
      } else {
        x = x_tilde;
        v = v_tilde;
      }
      if (update_modules) {
        module_update->updateModules(*this, rng);
      }
    }
#ifdef TIME_KEEPING
    void initialize_timers() {
      H_duration = std::chrono::duration<double>::zero();
      solver->DU_duration = std::chrono::duration<double>::zero();
      solver->approxDK_duration = std::chrono::duration<double>::zero();
    }
    template <typename StreamType>
    void print_timing_information(StreamType &stream) {
      stream << "---Sampling Timing Information" << std::endl;
      double DU_time = solver->DU_duration.count();
      double DK_time = solver->approxDK_duration.count();
      double H_time = H_duration.count();
      double total_time = H_time + DK_time + DU_time;
      stream << "Computing the Hamiltonian in time, " << H_time << " secs\n";
      stream << "Computing DU partial derivatives in time, " << DU_time
             << " secs\n";
      stream << "Computing DK partial derivatives in time, " << DK_time
             << " secs\n";
      stream << "H_time + DK_time + DU_time: " << total_time << "\n";
    }
#endif
  };
};

#endif
