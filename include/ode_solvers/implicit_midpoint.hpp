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
#include <chrono>

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
template <typename Point, typename NT, typename Polytope, typename func,int simdLen>
struct ImplicitMidpointODESolver {
  using VT = typename Polytope::VT;
  using MT = typename Polytope::MT;
  using pts = std::vector<MT>;
  using hamiltonian = Hamiltonian<Polytope, Point,simdLen>;
  using Opts = opts<NT>;

  unsigned int dim;
  int k=simdLen;
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

  hamiltonian ham;
  int num_runs=0;
  bool done;
#ifdef TIME_KEEPING
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> DU_duration =
      std::chrono::duration<double>::zero();
  std::chrono::duration<double> approxDK_duration =
      std::chrono::duration<double>::zero();
#endif
  ImplicitMidpointODESolver(NT initial_time, NT step, pts initial_state,
                            func oracle, Polytope &boundaries,
                            Opts &user_options)
      : eta(step), t(initial_time), xs(initial_state), F(oracle), P(boundaries),
        options(user_options), ham(hamiltonian(boundaries)) {
    dim = xs[0].rows();
  };

  void step(int k, bool accepted) {
    num_runs++;
    pts partialDerivatives;
#ifdef TIME_KEEPING
    start = std::chrono::system_clock::now();
#endif
if(num_runs<10){
std::cerr<<"---states------------\n";
std::cerr<< xs[0]<<"\n";
std::cerr<<"----------------------------\n";
std::cerr<< xs[1]<<"\n";
std::cerr<<"------------end------------\n";
}
partialDerivatives = ham.DU(xs);
    if(num_runs<10){
    std::cerr<<"partialDerivatives\n";
    std::cerr<< partialDerivatives[0]<<"\n";
    std::cerr<<"----------------------------\n";
    std::cerr<< partialDerivatives[1]<<"\n";
    std::cerr<<"---------end----------\n";
  }
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    DU_duration += end - start;
#endif
    xs = xs + partialDerivatives * (eta / 2);
    xs_prev = xs;
    done = false;
    nu = MT::Zero(P.equations(),simdLen);
    for (int i = 0; i < options.maxODEStep; i++) {
      pts xs_old = xs;
      pts xmid = (xs_prev + xs) / 2.0;
#ifdef TIME_KEEPING
      start = std::chrono::system_clock::now();
#endif
      // partialDerivatives = ham.DK(xmid);
      partialDerivatives = ham.approxDK(xmid, nu);
      if(num_runs<10){
      std::cerr<<"DK_partial\n";
      std::cerr<< partialDerivatives[0]<<"\n";
      std::cerr<<"----------------------------\n";
      std::cerr<< partialDerivatives[1]<<"\n";
      std::cerr<<"---------end----------\n";
    }
#ifdef TIME_KEEPING
      end = std::chrono::system_clock::now();
      approxDK_duration += end - start;
#endif
      xs = xs_prev + partialDerivatives * (eta);
      VT dist = ham.x_norm(xmid, xs - xs_old) / eta;
      NT maxdist = dist.maxCoeff();
      if(num_runs<10){
      std::cerr<<"===============states=================\n";
      std::cerr<< xs[0]<<"\n";
      std::cerr<<"----------------------------\n";
      std::cerr<< xs[1]<<"\n";
      std::cerr<<"------------end------------\n";
      std::cerr<<"maxdist= "<<maxdist<<"\n";
      }
      if (maxdist < options.implicitTol) {
        done = true;
        num_steps = i;
        break;
      } else if (maxdist > 1e16) {
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
    if(num_runs<10){
    std::cerr<<"---stateslast------------\n";
    std::cerr<< xs[0]<<"\n";
    std::cerr<<"----------------------------\n";
    std::cerr<< xs[1]<<"\n";
    std::cerr<<"------------end------------\n";
    }
    ham.project(xs);
  }

  void steps(int num_steps, bool accepted) {
    for (int i = 0; i < num_steps; i++) {
      step(i, accepted);
    }
  }

  MT get_state(int index) { return xs[index]; }

  void set_state(int index, MT p) { xs[index] = p; }
  void print_state() {
    for (int j = 0; j < xs.size(); j++) {
      std::cerr << "state " << j << ": ";
      for (unsigned int i = 0; i < xs[j].cols(); i++) {
        std::cerr << xs[j][i] << " ";
      }
      std::cerr << '\n';
    }
  }
};

#endif
