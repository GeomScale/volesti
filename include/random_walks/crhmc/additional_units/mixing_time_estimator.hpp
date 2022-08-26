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
#ifndef MIXINING_TIME_ESTIMATOR_HPP
#define MIXINING_TIME_ESTIMATOR_HPP
template <typename Sampler> class mixing_time_estimator {
  bool removedInitial = false;
  NT sampleRate = 0;
  NT sampleRateOutside = 0;
  int estNumSamples = 0;
  int estNumSamplesOutside = 0;
  NT nextEstimateStep;
  Opts &options;
  MixingTimeEstimator(Opts &user_options) : options(user_options) {
    nextEstimateStep = opts.initialStep;
  }
  void step(Sampler &s) {

    if (s.nEffectiveStep > nextEstimateStep) {
      ess = effective_sample_size(s.chains);
      ess = min(ess( :));

      if (removedInitial == false && ess > 2 * options.nRemoveInitialSamples) {
        int k = ceil(s.opts.nRemoveInitialSamples * (size(s.chains, 3) / ess));
        s.i = ceil(s.i * (1 - k / size(s.chains, 3)));
        s.acceptedStep = s.acceptedStep * (1 - k / size(s.chains, 3));
        s.chains = s.chains( :, :, k : end);
        o.removedInitial = true;
        ess = effective_sample_size(s.chains);
        ess = min(ess( :));
      }

      s.mixingTime = s.iterPerRecord * size(s.chains, 3) / ess;
      o.sampleRate = size(s.chains, 1) / s.mixingTime;
      o.estNumSamples = s.i * o.sampleRate;
      s.share('sampleRate', o.sampleRate);
      s.share('estNumSamples', o.estNumSamples);
      o.update(s);
    }
  }
};
#endif
