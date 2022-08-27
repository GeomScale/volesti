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
#include "Eigen/Eigen"
// Module for updating the extra term we add to the barrier
// This is nessecary for any polytope with free variables
template <typename Sampler, typename RandomNumberGenerator>
class dynamic_regularizer {
public:
  using NT = typename Sampler::NT;
  using Point = typename Sampler::point;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using Opts = typename Sampler::Opts;
  int n;
  VT bound;
  Opts &options;
  VT &extraHessian;
  dynamic_regularizer(Sampler &s)
      : options(s.params.options),
        extraHessian(options.DynamicWeight
                         ? s.solver->ham.weighted_barrier->extraHessian
                         : s.solver->ham.barrier->extraHessian) {
    n = s.dim;
    bound = VT::Ones(n);
    extraHessian = VT::Ones(n);
  }

  void update_regularization_factor(Sampler &s, RandomNumberGenerator &rng) {
    VT x = s.x.getCoefficients();
    x = (x.cwiseAbs()).cwiseMax(VT::Ones(n));
    bound = bound.cwiseMax(x);
    bool Condition =
        (2 / (bound.array() * bound.array()) < n * extraHessian.array()).any();

    if (Condition) {
      extraHessian = (0.5 / n) * (bound.cwiseProduct(bound)).cwiseInverse();
      s.solver->ham.move({s.x, s.v});
      s.v = s.GetDirectionWithMomentum(n, rng, s.x, Point(n), false);
    }
    /*
    for (int i = 0; i < n; i++) {
      bool change = false;
      if (options.DynamicWeight) {
        change = (2 / (bound(i) * bound(i))) <
                 n * extraHessian(i);
        if (change) {
          extraHessian =
              (0.5 / n) * (bound.cwiseProduct(bound)).cwiseInverse();
          s.solver->ham.move({s.x, s.v});
          s.v = s.GetDirectionWithMomentum(n, rng, s.x, Point(n), false);
          break;
        }
      } else {
        change = (2 / (bound(i) * bound(i))) <
                 n * s.solver->ham.barrier->extraHessian(i);
        if (change) {
          s.solver->ham.barrier->extraHessian =
              (0.5 / n) * (bound.cwiseProduct(bound)).cwiseInverse();
          s.solver->ham.move({s.x, s.v});
          s.v = s.GetDirectionWithMomentum(n, rng, s.x, Point(n), false);
          break;
        }
      }
    }
    */
  }
};
#endif
