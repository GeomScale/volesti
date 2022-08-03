// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef IMPLICIT_MIDPOINT_HPP
#define IMPLICIT_MIDPOINT_HPP

template <typename Point, typename NT, typename Polytope, typename func>
struct ImplicitMidpointODESolver {

  typedef std::vector<Point> pts;
  typedef std::vector<Polytope *> bounds;
  typedef typename Polytope::VT VT;
  typedef Hamiltonian<Point,Polytope> Ham;

  unsigned int dim;

  NT eta;
  NT t;

  func F;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  // Function oracle
  Hamiltonian<Point, Polytope> *ham;

  bool done;

  ImplicitMidpointODESolver(NT initial_time, NT step, pts initial_state,
                            func oracle, bounds boundaries, Ham *hamiltonian)
      : eta(step), t(initial_time), F(oracle),xs(initial_state), ham(hamiltonian) {
    dim = xs[0].dimension();
  };
  pts operator+(const pts &v1, const pts &v2) {
    pts result;
    for (int i = 0; i < v1.size(); i++) {
      result[i] = v1[i] + v2[i];
    }
    return result;
  }
  pts operator*(const pts &v1, NT alfa) {
    pts result;
    for (int i = 0; i < v1.size(); i++) {
      result[i] = v1[i] * alfa;
    }
    return result;
  }
  pts operator/(const pts &v1, NT alfa) { return v1 * (1 / alfa); }

  void step(int k, bool accepted) {
    pts xs_old = xs;
    xmid = (xs_prev + xs) / 2;
    pts partialDerivatives = ham->approxDK(xmid, vmid);
    xs = xs_prev + eta * partialDerivatives;
    VT dist = ham->x_norm(xmid[0], xs[0] - xs_old[0]) / eta;
    NT maxdist = dist.maxCoeff();
    if (maxdist < implicitTol) {
      done = true;
    }
  }

  void steps(int num_steps, bool accepted) {
    xs = xs + (-h / 2) * ham->DU(xs[0]);
    xs_prev = xs;
    done = false;
    for (int i = 0; i < num_steps; i++) {
      if (done) {
        break;
      }
      step(i, accepted);
    }
    xs = xs + (-h / 2) * ham->DU(xs[0]);
    ham->prepare(xs[0]);
    xs = ham->project(xs[0]);
  }

  Point get_state(int index) { return xs[index]; }

  void set_state(int index, Point p) { xs[index] = p; }
};

#endif
