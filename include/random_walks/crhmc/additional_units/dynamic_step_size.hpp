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
template <typename Sampler> class dynamic_step_size {
  int consecutiveBadStep = 0;
  int iterSinceShrink = 0;
  NT rejectSinceShrink = 0;
  int ODEStepSinceShrink = 0;
  int effectiveStep = 0;
  warmupFinished = false;
  Opts &options;
  dynamic_step_size(Sampler &s) : options(s.params.options) {
    if (warmUpStep > 0) {
      s.solver->eta = 1e-3;
    } else {
      warmupFinished = true;
    }
  }
  udate_step_size(Sampler &s) {
    int bad_step = 0;
    eta = s->solver.eta;
    if (s.prob < 0.5 || s.ode_step == s.options.maxODEStep) {
      bad_step = 1;
    }
    consecutiveBadStep = bad_step * consecutiveBadStep + bad_step;
    NT warmupRatio = s.auto_tuner.nEffectiveStep / options.warmUpStep;
    if (warmupRatio < 1 && !warmupFinished &&
        consecutiveBadStep < options.maxConsecutiveBadStep) {
      s->solver.eta = options.initalStepSize * std::min(warmupRatio + 1e-2, 1);
      s.params.momentum =
          1 - std::min(1.0, s->solver.eta / options.effectiveStepSize);
      return;
    }
    if (!warmupFinished) {
      s.i = 1;
      s.acceptedStep = 0;
      s.nEffectiveStep = 0;
      warmupFinished = true;
    }

    iterSinceShrink = iterSinceShrink + 1;
    rejectSinceShrink = rejectSinceShrink + 1 - s.prob;
    ODEStepSinceShrink = ODEStepSinceShrink + s.ode_step;

    shrink = 0;
    shiftedIter = iterSinceShrink + 20 / (1 - s.params.momentum);

    targetProbability = (1 - s.params.momentum) ^ (2 / 3) / 4;
    if (rejectSinceShrink > targetProbability * shiftedIter) {
      shrink = 1;
    }

    if (consecutiveBadStep > options.maxConsecutiveBadStep) {
      shrink = 1;
    }

    if (ODEStepSinceShrink > options.targetODEStep * shiftedIter) {
      shrink = 1;
    }

    if (shrink) {
      iterSinceShrink = 0;
      rejectSinceShrink = 0;
      ODEStepSinceShrink = 0;
      consecutiveBadStep = 0;

      s->solver.eta /= options.shrinkFactor;
      s.params.momentum =
          1 - std::min(0.999, s->solver.eta / options.effectiveStepSize);

      if (s->solver.eta < options.minStepSize) {
        std::cerr << "Algorithm fails to converge even with step size h = "
                  << eta << "\n";
        exit(1);
      }
    }

    iterSinceShrink = iterSinceShrink + 1;
  }
};
#endif
