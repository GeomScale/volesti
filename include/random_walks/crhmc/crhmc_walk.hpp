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
    using IVT = Eigen::Matrix<int, Eigen::Dynamic, 1>;

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
    auto_tuner<Sampler, RandomNumberGenerator> *module_update;

    // Helper variables
    VT H, H_tilde;
    // Density exponent
    NegativeLogprobFunctor &f;
#ifdef TIME_KEEPING
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> H_duration =
        std::chrono::duration<double>::zero();
#endif
    Walk(Polytope &Problem, MT &p, NegativeGradientFunctor &neg_grad_f,
         NegativeLogprobFunctor &neg_logprob_f,
         parameters<NT, NegativeGradientFunctor> &param)
        : params(param), F(neg_grad_f), f(neg_logprob_f), P(Problem) {

      dim = p.rows();

      // Starting point is provided from outside
      x = p;
      accepted = false;
      // Initialize solver
      solver =
          new Solver(0.0, params.eta, {x, x}, F, Problem, params.options);
      simdLen=solver->k;
      v = solver->get_state(1);
      module_update = new auto_tuner<Sampler, RandomNumberGenerator>(*this);
      update_modules = params.options.DynamicWeight ||
                       params.options.DynamicRegularizer ||
                       params.options.DynamicStepSize;
    };
    //Sample a new velocity with momentum
    MT GetDirectionWithMomentum(unsigned int const &dim,
                                   RandomNumberGenerator &rng, MT x, MT v,
                                   NT momentum = 0, bool normalize = true) {
      MT z=MT(dim,simdLen);
      for(int i=0;i<simdLen;i++){
        z(Eigen::all,i)=GetDirection<Point>::apply(dim, rng, normalize).getCoefficients();
      }
      solver->ham.move({x, v});
      MT sqrthess = (solver->ham.hess).cwiseSqrt();
      z = sqrthess.cwiseProduct(z);
      return v * std::sqrt(momentum) + z * std::sqrt(1 - momentum);
    }
    // Returns the current point in the tranformed in the original space
    inline MT getPoints() {
      return (P.T * x).colwise() + P.y;
     }

    // Returns the current point in the tranformed in the original space
    inline Point getPoint() { return Point(P.T * x + P.y); }
    inline MT masked_choose(MT &x,MT &x_tilde,IVT &accept){
      MT result=MT(x.rows(),x.cols());
      for(int i=0;i<simdLen;i++){
        if(accept(i)==1){
          result(Eigen::all,i)=x_tilde(Eigen::all,i);
        }else{
          result(Eigen::all,i)=x(Eigen::all,i);
        }
      }
      return result;
    }
    inline void apply(RandomNumberGenerator &rng, int walk_length = 1,
                      bool metropolis_filter = true) {

      num_runs++;
      //  Pick a random velocity with momentum
      v = GetDirectionWithMomentum(dim, rng, x, v, params.momentum, false);
      if(num_runs<10){
        std::cerr<<"x=\n"<<x<<"\n";
        std::cerr<<"v=\n"<<v<<"\n";
      }
      solver->set_state(0, x);
      solver->set_state(1, v);
      // Get proposals
      solver->steps(walk_length, accepted);
      x_tilde = solver->get_state(0);
      v_tilde = solver->get_state(1);
      if(num_runs<10){
        std::cerr<<"x_tilde=\n"<<x_tilde<<"\n";
        std::cerr<<"v_tilde=\n"<<v_tilde<<"\n";
      }

      if (metropolis_filter) {
#ifdef TIME_KEEPING
        start = std::chrono::system_clock::now();
#endif
        // Calculate initial Hamiltonian
        H = solver->ham.hamiltonian(x, v);

        // Calculate new Hamiltonian
        H_tilde = solver->ham.hamiltonian(x_tilde, - v_tilde);

#ifdef TIME_KEEPING
        end = std::chrono::system_clock::now();
        H_duration += end - start;
#endif
        VT feasible = solver->ham.feasible(x_tilde,
                                           v_tilde);
        prob=(1.0< exp((H - H_tilde).array())).select(1.0,exp((H - H_tilde).array()));

        prob=prob.cwiseProduct(feasible);
        if(num_runs<10){
        std::cerr<<"--------prob*feasible----------\n"<<prob<<"\n";
        }
        total_acceptance_prob += prob.sum();
        VT rng_vector=VT(simdLen);
        for(int i=0;i<simdLen;i++){rng_vector(i)=rng.sample_urdist();}
        accept=(rng_vector.array()<prob.array()).select(1*IVT::Ones(simdLen),0*IVT::Ones(simdLen));
        
        x=masked_choose(x,x_tilde,accept);
        v=-v;
        v=masked_choose(v,v_tilde,accept);
        discard_ratio = (1.0 * total_discarded_samples) / num_runs;
        average_acceptance_prob = total_acceptance_prob / num_runs;
      } else {
        x = x_tilde;
        v = v_tilde;
      }
      if(num_runs<10){
        std::cerr<<"x_new=\n"<<x<<"\n";
        std::cerr<<"v_new=\n"<<v<<"\n";
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
      std::cerr << "Computing the Hamiltonian in time, " << H_time << " secs\n";
      std::cerr << "Computing DU partial derivatives in time, " << DU_time
                << " secs\n";
      std::cerr << "Computing DK partial derivatives in time, " << DK_time
                << " secs\n";
      std::cerr << "H_time + DK_time + DU_time: " << total_time << "\n";
    }
#endif
  };
};

#endif
