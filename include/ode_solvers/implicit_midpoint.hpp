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

#ifndef IMPLICIT_MIDPOINT_HPP
#define IMPLICIT_MIDPOINT_HPP
#include "preprocess/crhmc/opts.h"
#include "random_walks/crhmc/hamiltonian.hpp"

template <typename T>
inline std::vector<T> operator+(const std::vector<T> &v1,
                                const std::vector<T> &v2) {
  std::vector<T> result(v1.size());
  for (int i = 0; i < v1.size(); i++) {
    result[i] = v1[i] + v2[i];
  }
  return result;
}
template <typename T, typename Type>
inline std::vector<T> operator*(const std::vector<T> &v, const Type alfa) {
  std::vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = v[i] * alfa;
  }
  return result;
}
template <typename T, typename Type>
inline std::vector<T> operator/(const std::vector<T> &v, const Type alfa) {
  return v * (1 / alfa);
}
template <typename T>
inline std::vector<T> operator-(const std::vector<T> &v1,
                                const std::vector<T> &v2) {

  return v1 + v2 * (-1.0);
}
template <typename Point, typename NT, typename Polytope, typename func>
struct ImplicitMidpointODESolver {
  using VT = typename Polytope::VT;
  using MT = typename Polytope::MT;
  using pts = std::vector<Point>;
  using hamiltonian = Hamiltonian<Polytope, Point, func>;
  using Opts = opts<NT>;

  unsigned int dim;

  NT eta;
  NT t;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  // Function oracle
  func F;
  Polytope &P;
  Opts& options;
  MT nu;

  hamiltonian ham;

  bool done;

  ImplicitMidpointODESolver(NT initial_time, NT step, pts initial_state,
                            func oracle, Polytope &boundaries,
                            Opts &user_options)
      : eta(step), t(initial_time), xs(initial_state), F(oracle),
        options(user_options), P(boundaries),
        ham(hamiltonian(boundaries, oracle)) {
    dim = xs[0].dimension();
  };

  void step(int k, bool accepted) {
    pts xs_old = xs;
    pts xmid = (xs_prev + xs) / 2.0;
    pts partialDerivatives;
    partialDerivatives = ham.DK(xmid);
    xs = xs_prev + partialDerivatives * (eta);
    VT dist = ham.x_norm(xmid, xs - xs_old) / eta;
    NT maxdist = dist.maxCoeff();
    if (maxdist < options.implicitTol) {
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
    xs[0] = Point(P.project(xs[0].getCoefficients()));
  }

  Point get_state(int index) { return xs[index]; }

  void set_state(int index, Point p) { xs[index] = p; }
  void print_state() {
    for (int j = 0; j < xs.size(); j++) {
      std::cout << "state " << j << ": ";
      for (unsigned int i = 0; i < xs[j].dimension(); i++) {
        std::cout << xs[j][i] << " ";
      }
      std::cout << '\n';
    }
  }
};

#endif
