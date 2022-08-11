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
template <typename NT> class mixing_time_estimator {
  using Opts = opts<NT>;
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
  void step() {}
};
#endif
