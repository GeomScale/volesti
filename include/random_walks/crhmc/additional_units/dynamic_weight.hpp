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

/// Class responsible for updating the weights of the barrier
/// Part of crhmc sampler
template <typename Sampler, typename RandomNumberGenerator>
class dynamic_weight {
  using NT = typename Sampler::NT;
  using Point = typename Sampler::point;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using IVT = Eigen::Array<int, Eigen::Dynamic, 1>;
  using Opts = typename Sampler::Opts;

public:
  int simdLen;
  IVT consecutiveBadStep;
  int n;
  VT &w;
  Opts options;
  dynamic_weight(Sampler &s)
      : simdLen(s.simdLen), w(s.solver->ham.weighted_barrier->w), options(s.params.options)
  {
    n = s.dim;
    consecutiveBadStep = IVT::Zero(simdLen);
  }
  // If we have consecutive bad steps update the weights with
  //  the help of the leverage scores.
  void update_weights(Sampler &s, RandomNumberGenerator &rng)
  {
    IVT bad_step = IVT::Zero(simdLen);
    if (s.solver->num_steps == options.maxODEStep) {
      bad_step += 1;
    } else {
      bad_step = (s.prob.array() < 0.5).select(1, IVT::Zero(simdLen));
    }
    NT threshold;
    consecutiveBadStep = bad_step * consecutiveBadStep + bad_step;

    if (s.accept.sum() < simdLen) {
      VT lsc = s.solver->ham.lsc.colwise().maxCoeff().transpose();
      /*The more bad steps in a row we have the higher we want the threshold to be
      In order to change w more drasticaly according to the leverage scores.
      So if we have more than 2 bad steps in a row we elect to set the threshold to 4
      else to 16. Not many changes will be possible as the w should be upperbounded by 1*/
      if (consecutiveBadStep.maxCoeff() > 2) {
        threshold = 4;
      } else {
        threshold = 16;
      }
      bool changed = (lsc.array() > threshold * w.array()).any();
      if (changed) {
        w = (lsc.array() > threshold * w.array()).select((w * threshold).cwiseMin(1), w);
        s.solver->ham.forceUpdate = true;
        s.solver->ham.move({s.x, s.v});
        s.v = s.get_direction_with_momentum(n, rng, s.x, MT::Zero(n, simdLen), false);
      }
    }
  }
};
#endif
