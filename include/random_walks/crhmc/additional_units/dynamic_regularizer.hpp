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

/// Module for updating the extra term we add to the barrier
/// This is nessecary for any polytope with free variables
/// Part of crhmc sampler
template <typename Sampler, typename RandomNumberGenerator>
class dynamic_regularizer {
public:
  using NT = typename Sampler::NT;
  using Point = typename Sampler::point;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using Opts = typename Sampler::Opts;
  int n;
  int simdLen;
  MT bound;
  Opts &options;
  MT &extraHessian;
  dynamic_regularizer(Sampler &s) :
    simdLen(s.simdLen),
    options(s.params.options),
    extraHessian(options.DynamicWeight
      ? s.solver->ham.weighted_barrier->extraHessian
      : s.solver->ham.barrier->extraHessian)
  {
    n = s.dim;
    bound = MT::Ones(n, simdLen);
    extraHessian = MT::Ones(n, simdLen);
  }

  void update_regularization_factor(Sampler &s, RandomNumberGenerator &rng) {
    MT x = s.x;
    x = (x.cwiseAbs()).cwiseMax(1);
    bound = bound.cwiseMax(x);
    if ((2 / (bound.array() * bound.array()) < n * extraHessian.array()).any()) {
      extraHessian = (0.5 / n) * (bound.cwiseProduct(bound)).cwiseInverse();
      s.solver->ham.forceUpdate = true;
      s.solver->ham.move({s.x, s.v});
      s.v = s.get_direction_with_momentum(n, rng, s.x, MT::Zero(n, simdLen), 0, false);
    }
  }
};
#endif
