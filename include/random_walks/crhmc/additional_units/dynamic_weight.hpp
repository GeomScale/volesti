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
#ifndef DYNAMIC_WEIGHT_HPP
#define DYNAMIC_WEIGHT_HPP
template <typename Sampler> class dynamic_weight {
public:
  int consecutiveBadStep = 0;
  DynamicWeight(Sampler &s) {}
  step(Sampler &s) {
    int bad_step = 0;
    if (s.prob < 0.5 || s.ode_step == s.options.maxODEStep) {
      bad_step = 1;
    }
    consecutiveBadStep = bad_step * consecutiveBadStep + bad_step;

    if (!s.freezed && !s.accept) {
      VT lsc =s->solver.ham.lsc;
      if (consecutiveBadStep > 2) {
        threshold = 4;
      } else {
        threshold = 16;
      }
      int n = w.rows();
      bool changed = false;
      for (int i = 0; i < n; i++) {
        if (lsc(i) > threshold * s->solver.ham->barrier.w(i)) {
          s.ham.barrier.w(i) = std::min(s.ham.barrier.w(i) * threshold, 1);
          changed = 1;
        }
      }
      if(changed)
        s.ham.move(s.solver.xs, true);
    }
  }
};
#endif
