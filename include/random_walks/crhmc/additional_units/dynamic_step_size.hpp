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

public:
  int consecutiveBadStep = 0;
  int iterSinceShrink = 0;
  NT rejectSinceShrink = 0;
  int ODEStepSinceShrink = 0;
  int effectiveStep = 0;
  bool warmupFinished = false;
  Opts &options;
  NT &eta;
  NT &momentum;
  NT acceptedStep = 0;
  NT accumulatedMomentum = 0;
  NT nEffectiveStep = 0; // number of effective steps

  dynamic_step_size(Sampler &s) :
    options(s.params.options),
    eta(s.solver->eta),
    momentum(s.params.momentum)
  {
    if (options.warmUpStep > 0) {
      eta = 1e-3;
    } else {
      warmupFinished = true;
    }
  }
  void update_step_size(Sampler &s) {
    acceptedStep = acceptedStep + s.prob;
    accumulatedMomentum = s.prob * momentum * accumulatedMomentum + eta;
    nEffectiveStep = nEffectiveStep + eta * accumulatedMomentum * s.accept;

    int bad_step = s.prob < 0.5 || s.solver->num_steps == options.maxODEStep ? 1 : 0;
    consecutiveBadStep = bad_step * consecutiveBadStep + bad_step;

    NT warmupRatio = nEffectiveStep / options.warmUpStep;
    if (warmupRatio < 1 && !warmupFinished &&
        consecutiveBadStep < options.maxConsecutiveBadStep) {
      eta = options.initialStep * std::min(warmupRatio + 1e-2, 1.0);
      momentum = 1 - std::min(1.0, eta / options.effectiveStepSize);
      return;
    }
    if (!warmupFinished) {
      acceptedStep = 0;
      nEffectiveStep = 0;
      warmupFinished = true;
    }

    iterSinceShrink++;
    rejectSinceShrink += 1 - s.prob;
    ODEStepSinceShrink += s.solver->num_steps;

    int shrink = 0;
    NT shiftedIter = iterSinceShrink + 20 / (1 - momentum);

    NT targetProbability = std::pow((1.0 - momentum), (2 / 3)) / 4;
    if (rejectSinceShrink > targetProbability * shiftedIter ||
        consecutiveBadStep > options.maxConsecutiveBadStep  ||
        ODEStepSinceShrink > options.targetODEStep * shiftedIter) {
      shrink = 1;
    }

    if (shrink == 1) {
      iterSinceShrink = 0;
      rejectSinceShrink = 0;
      ODEStepSinceShrink = 0;
      consecutiveBadStep = 0;

      eta /= options.shrinkFactor;
      momentum = 1 - std::min(0.999, eta / options.effectiveStepSize);

      if (eta < options.minStepSize) {
        std::cerr << "Algorithm fails to converge even with step size h = "
                  << eta << "\n";
        exit(1);
      }
    }

    iterSinceShrink = iterSinceShrink + 1;
  }
};
#endif
