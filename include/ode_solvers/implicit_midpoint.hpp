// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef IMPLICIT_MIDPOINT_HPP
#define IMPLICIT_MIDPOINT_HPP
#include "preprocess/crhmc/opts.h"
#include "random_walks/crhmc/hamiltonian_utils.hpp"

template <typename T>
std::vector<T> operator+(const std::vector<T> &v1, const std::vector<T> &v2) {
  std::vector<T> result(v1.size());
  for (int i = 0; i < v1.size(); i++) {
    result[i] = v1[i] + v2[i];
  }
  return result;
}
template <typename T, typename Type>
std::vector<T> operator*(const std::vector<T> &v, const Type alfa) {
  std::vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = v[i] * alfa;
  }
  return result;
}
template <typename T, typename Type>
std::vector<T> operator/(const std::vector<T> &v, const Type alfa) {
  return v * (1 / alfa);
}

template <typename Point, typename NT, typename Polytope, typename func>
struct ImplicitMidpointODESolver {

  typedef typename Polytope::VT VT;
  typedef typename Polytope::VT MT;
  typedef std::vector<VT> pts;
  typedef Hamiltonian<Polytope, func> hamiltonian;
  using Opts = opts<NT>;

  unsigned int dim;

  NT eta;
  NT t;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  // Function oracle
  func F;
  Polytope *P;
  Opts *options;
  MT nu;

  hamiltonian ham;

  bool done;

  ImplicitMidpointODESolver(NT initial_time, NT step, pts initial_state,
                            func oracle, Polytope *boundaries,
                            Opts *user_options)
      : eta(step), t(initial_time), xs(initial_state), F(oracle),
        options(user_options), P(boundaries),
        ham(hamiltonian(boundaries, oracle)) {
    dim = xs[0].rows();
  };

  void step(int k, bool accepted) {
    pts xs_old = xs;
    pts xmid = (xs_prev + xs) / 2.0;
    pts partialDerivatives;
    partialDerivatives = ham.DK(xmid);
    xs = xs_prev + partialDerivatives * (eta);
    VT dist = ham.x_norm(xmid[0], xs[0] - xs_old[0]) / eta;
    NT maxdist = dist.maxCoeff();
    if (maxdist < options->implicitTol) {
      done = true;
    }
  }

  void steps(int num_steps, bool accepted) {
    pts partialDerivatives = ham.DU(xs);
    xs = xs + partialDerivatives * (eta / 2);
    xs_prev = xs;
    done = false;
    for (int i = 0; i < num_steps; i++) {
      if (done) {
        break;
      }
      step(i, accepted);
    }
    partialDerivatives = ham.DU(xs);
    xs = xs + partialDerivatives * (eta / 2);
    xs[0] = P->project(xs[0]);
  }

  Point get_state(int index) { return xs[index]; }

  void set_state(int index, Point p) { xs[index] = p; }
  void print_state() {
    for (int j = 0; j < xs.size(); j++) {
      std::cout << xs[j].transpose() << std::endl;
    }
    std::cout << std::endl;
  }
};

#endif
