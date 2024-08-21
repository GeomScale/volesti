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
#ifdef TIME_KEEPING
#include <chrono>
#endif

template <typename T>
inline std::vector<T> operator+(const std::vector<T> &v1,
                                const std::vector<T> &v2)
{
  std::vector<T> result(v1.size());
  for (int i = 0; i < v1.size(); i++) {
    result[i] = v1[i] + v2[i];
  }
  return result;
}
template <typename T, typename Type>
inline std::vector<T> operator*(const std::vector<T> &v, const Type alpha)
{
  std::vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = v[i] * alpha;
  }
  return result;
}
template <typename T, typename Type>
inline std::vector<T> operator/(const std::vector<T> &v, const Type alpha)
{
  return v * (1 / alpha);
}
template <typename T>
inline std::vector<T> operator-(const std::vector<T> &v1,
                                const std::vector<T> &v2)
{

  return v1 + v2 * (-1.0);
}
template <typename Point, typename NT, typename Polytope, typename func, int simdLen = 1>
struct ImplicitMidpointODESolver {
  using VT = typename Polytope::VT;
  using MT = typename Polytope::MT;
  using pts = std::vector<MT>;
  using hamiltonian = Hamiltonian<Polytope, Point, simdLen>;
  using Opts = opts<NT>;

  unsigned int dim;
  NT eta;
  int num_steps = 0;
  NT t;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  // Function oracle
  func F;
  Polytope &P;
  Opts &options;
  MT nu;
  int num_runs = 0;
  hamiltonian ham;
  bool done;
#ifdef TIME_KEEPING
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> DU_duration = std::chrono::duration<double>::zero();
  std::chrono::duration<double> approxDK_duration = std::chrono::duration<double>::zero();
#endif
  ImplicitMidpointODESolver(NT initial_time,
                            NT step,
                            pts initial_state,
                            func oracle,
                            Polytope &boundaries,
                            Opts &user_options) :
                            eta(step),
                            t(initial_time),
                            xs(initial_state),
                            F(oracle),
                            P(boundaries),
                            options(user_options),
                            ham(hamiltonian(boundaries))
  {
    dim = xs[0].rows();
  };

  void step(int k, bool accepted) {
    num_runs++;
    pts partialDerivatives;
#ifdef TIME_KEEPING
    start = std::chrono::system_clock::now();
#endif
    partialDerivatives = ham.DU(xs);
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    DU_duration += end - start;
#endif
    xs = xs + partialDerivatives * (eta / 2);
    xs_prev = xs;
    done = false;
    nu = MT::Zero(P.equations(), simdLen);
    for (int i = 0; i < options.maxODEStep; i++) {
      pts xs_old = xs;
      pts xmid = (xs_prev + xs) / 2.0;
#ifdef TIME_KEEPING
      start = std::chrono::system_clock::now();
#endif
      partialDerivatives = ham.approxDK(xmid, nu);

#ifdef TIME_KEEPING
      end = std::chrono::system_clock::now();
      approxDK_duration += end - start;
#endif
      xs = xs_prev + partialDerivatives * (eta);
      VT dist = ham.x_norm(xmid, xs - xs_old) / eta;
      NT maxdist = dist.maxCoeff();
      //If the estimate does not change terminate
      if (maxdist < options.implicitTol) {
        done = true;
        num_steps = i;
        break;
      //If the estimate is very bad, sample another velocity
    } else if (maxdist > options.convergence_bound) {
        xs = xs * std::nan("1");
        done = true;
        num_steps = i;
        break;
      }
    }
#ifdef TIME_KEEPING
    start = std::chrono::system_clock::now();
#endif
    partialDerivatives = ham.DU(xs);
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    DU_duration += end - start;
#endif
    xs = xs + partialDerivatives * (eta / 2);
    ham.project(xs);
  }

  void steps(int num_steps, bool accepted) {
    for (int i = 0; i < num_steps; i++) {
      step(i, accepted);
    }
  }

  MT get_state(int index) { return xs[index]; }

  void set_state(int index, MT p) { xs[index] = p; }
  template<typename StreamType>
  void print_state(StreamType &stream) {
    for (int j = 0; j < xs.size(); j++) {
      stream << "state " << j << ": \n";
      stream << xs[j];
      stream << '\n';
    }
  }

  NT get_eta() const { return eta; }
};

#endif
