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
#ifndef DYNAMIC_REGULARIZER_HPP
#define DYNAMIC_REGULARIZER_HPP
// Module for updating the extra term we add to the barrier
// This is nessecary for any polytope with free variables
template <typename Sampler> void DynamicRegularizerUpdate(Sampler s) {
  VT x = s.x.getCoefficients();
  int n = x.rows();
  x = (x.cwiseAbs()).cwiseMax(VT::Ones(n));
  bound = bound.cwiseMax(x);
  if (!s.freezed) {
    bool changed = false;
    for (int i = 0; i < n; i++) {
      bool change = false;
      if (options.DynamicWeight) {
        change = (2 / (bound(i) * bound(i))) <
                 n * s->solver.ham->weighted_barrier.extraHessian(i);
        if (change) {
          s->solver.ham->weighted_barrier.extraHessian =
              0.5 / n * (bound.cwiseProduct(bound)).cwiseInverse();
          s->solver.ham.move(s->solver.xs, true);
        }

      } else {
        change = (2 / (bound(i) * bound(i))) <
                 n * s->solver.ham->barrier.extraHessian(i);
        if (change) {
          s->solver.ham->barrier.extraHessian =
              0.5 / n * (bound.cwiseProduct(bound)).cwiseInverse();
          s->solver.ham.move(s->solver.xs, true);
        }
      }
      if (change) {
        break;
      }
    }
  }
}
#endif
