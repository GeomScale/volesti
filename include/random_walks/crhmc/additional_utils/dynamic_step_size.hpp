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
template <typename NT> class dynamic_step_size {
  int consecutiveBadStep = 0;
  int iterSinceShrink = 0;
  int rejectSinceShrink = 0;
  int ODEStepSinceShrink = 0;
  int effectiveStep = 0;
  warmupFinished = false;
  dynamic_step_size(s) {
    if (warmUpStep > 0) {
      stepSize = 1e-3;
    } else {
      warmupFinished = true;
    }
  }
};
#endif
