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
#ifndef DYNAMIC_STEP_SIZE_HPP
#define DYNAMIC_STEP_SIZE_HPP

/// Module for dynamically choosing the ODE step size and the velocity momentum
/// Part of crhmc sampler
template <typename Sampler>
class dynamic_step_size {
  using NT = typename Sampler::NT;
  using Opts = typename Sampler::Opts;
  using IVT = Eigen::Array<int, Eigen::Dynamic, 1>;
  using VT = Eigen::Array<NT, Eigen::Dynamic, 1>;

public:
  int simdLen;
  IVT consecutiveBadStep;
  int iterSinceShrink = 0;
  VT rejectSinceShrink;
  int ODEStepSinceShrink = 0;
  int effectiveStep = 0;
  bool warmupFinished = false;
  Opts &options;
  NT &eta;
  NT &momentum;
  VT acceptedStep;
  VT nEffectiveStep; // number of effective steps
  NT accumulatedMomentum = 0;
  dynamic_step_size(Sampler &s) :
    simdLen(s.simdLen),
    options(s.params.options),
    eta(s.solver->eta),
    momentum(s.params.momentum)
  {
    nEffectiveStep = VT::Zero(simdLen);
    acceptedStep = VT::Zero(simdLen);
    consecutiveBadStep = IVT::Zero(simdLen);
    rejectSinceShrink = VT::Zero(simdLen);

    if (options.warmUpStep > 0 && options.etaInitialize) {
      eta = 1e-3;
    } else {
      warmupFinished = true;
    }
  }
  void update_step_size(Sampler &s) {
    acceptedStep = acceptedStep + s.prob.array();
    accumulatedMomentum = s.prob.mean() * momentum * accumulatedMomentum + eta;
    Eigen::Matrix<NT, Eigen::Dynamic, 1> accept = (s.accept.template cast<NT>());
    nEffectiveStep = nEffectiveStep + eta * accumulatedMomentum * accept.array();
    IVT bad_step = IVT::Zero(simdLen);
    if (s.solver->num_steps == options.maxODEStep) {
      bad_step += 1;
    } else {
      bad_step = (s.prob.array() < 0.5).select(1, IVT::Zero(simdLen));
    }
    consecutiveBadStep = bad_step * consecutiveBadStep + bad_step;

    NT warmupRatio = nEffectiveStep.mean() / options.warmUpStep;
    if (warmupRatio < 1 && !warmupFinished &&
        consecutiveBadStep.maxCoeff() < options.maxConsecutiveBadStep) {
      eta = options.initialStep * std::min(warmupRatio + 1e-2, 1.0);
      momentum = 1 - std::min(1.0, eta / options.effectiveStepSize);
      return;
    }

    if (!warmupFinished) {
      acceptedStep = VT::Zero(simdLen);
      nEffectiveStep = VT::Zero(simdLen);
      warmupFinished = true;
    }

    iterSinceShrink++;
    rejectSinceShrink += 1 - s.prob.array();
    ODEStepSinceShrink += s.solver->num_steps;

    int shrink = 0;
    NT shiftedIter = iterSinceShrink + 20 / (1 - momentum);

    NT targetProbability = std::pow((1.0 - momentum), (2 / 3)) / 4;
    if (rejectSinceShrink.maxCoeff() > targetProbability * shiftedIter||
        consecutiveBadStep.maxCoeff() > options.maxConsecutiveBadStep ||
        ODEStepSinceShrink > options.targetODEStep * shiftedIter) {
      shrink = 1;
    }

    if (shrink == 1) {
      iterSinceShrink = 0;
      ODEStepSinceShrink = 0;
      rejectSinceShrink = VT::Zero(simdLen);
      consecutiveBadStep = IVT::Zero(simdLen);

      eta /= options.shrinkFactor;
      momentum = 1 - std::min(0.999, eta / options.effectiveStepSize);

      if (eta < options.minStepSize) {
        s.P.terminate=true;
        s.P.terminate_message="Algorithm fails to converge even with step size h = "+std::to_string(eta)+"\n";
      }
    }

    iterSinceShrink = iterSinceShrink + 1;
  }
};
#endif
