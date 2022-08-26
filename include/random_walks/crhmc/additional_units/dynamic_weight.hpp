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
/*Class responsible for updating the weights of the barrier*/
template <typename Sampler, typename RandomNumberGenerator>
class dynamic_weight {
  using NT = typename Sampler::NT;
  using Point = typename Sampler::point;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

public:
  int consecutiveBadStep = 0;
  int n;
  VT &w;
  dynamic_weight(Sampler &s) : w(s.solver->ham.weighted_barrier->w) {
    std::cerr << "Using dynamic weights------------------------------------"
              << '\n';
    n = s.dim;
  }
  // If we have consecutive bad steps update the weights with
  //  the help of the leverage scores.
  void update_weights(Sampler &s, RandomNumberGenerator &rng) {
    int bad_step = 0;
    NT threshold;
    if (s.prob < 0.5 || s.solver->num_steps == s.params.options.maxODEStep) {
      bad_step = 1;
    }
    consecutiveBadStep = bad_step * consecutiveBadStep + bad_step;

    if (!s.accepted) {
      VT lsc = s.solver->ham.lsc;
      if (consecutiveBadStep > 2) {
        threshold = 4;
      } else {
        threshold = 16;
      }
      bool changed = (lsc.array() > threshold * w.array()).any();
      w = (lsc.array() > threshold * w.array())
              .select((w * threshold).cwiseMin(VT::Ones(n)), w);
      //      (w.array() < lsc.array())
      //          .select(
      //              std::min(w(i) *
      //              threshold, 1.0), w);
      /*
            for (int i = 0; i < n; i++) {
              if (lsc(i) > threshold * w(i)) {
                w(i) = std::min(w(i) * threshold, 1.0);
                changed = true;
              }
            }
            */
      if (changed) {
        s.solver->ham.forceUpdate = true;
        s.solver->ham.move({s.x, s.v});
        s.v = s.GetDirectionWithMomentum(n, rng, s.x, Point(n), false);
      }
    }
  }
};
#endif
