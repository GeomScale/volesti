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

// The log barrier for the domain {lu <= x <= ub}:
//	phi(x) = - sum log(x - lb) - sum log(ub - x).
#ifndef TWOSIDEDBARIER_H
#define TWOSIDEDBARIER_H

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include <vector>

template <typename Point> class TwoSidedBarrier {

  using NT = typename Point::FT;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

public:
  VT lb;
  VT ub;
  int vdim;
  int n;
  std::vector<int> upperIdx;
  std::vector<int> lowerIdx;
  std::vector<int> freeIdx;
  VT center;
  const NT max_step = 1e16; // largest step size
  VT extraHessian;

  const NT inf = std::numeric_limits<NT>::infinity();

  void set_bound(VT const &_lb, VT const &_ub) {

    lb = _lb;
    ub = _ub;
    n = lb.rows();
    extraHessian = (1e-20) * VT::Ones(n);
    int x1 = 0, x2 = 0, x3 = 0;
    for (int i = 0; i < n; i++) {
      if (lb(i) == -inf) {
        upperIdx.push_back(i);
        x1++;
      }
      if (ub(i) == inf) {
        lowerIdx.push_back(i);
        x2++;
      }
      if (ub(i) == inf && lb(i) == -inf) {
        freeIdx.push_back(i);
        x3++;
      }
    }

    VT c = (ub + lb) / 2;

    c(lowerIdx) = lb(lowerIdx) + VT::Ones(x2, 1) * 1e6;
    c(upperIdx) = ub(upperIdx) - VT::Ones(x1, 1) * 1e6;
    c(freeIdx) *= 0.0;

    center = c;
  }
  TwoSidedBarrier(VT const &_lb, VT const &_ub, int _vdim = 1) {
    set_bound(_lb, _ub);
    vdim = _vdim;
    extraHessian = (1e-20) * VT::Ones(n);
  }
  TwoSidedBarrier() { vdim = 1; }

  VT gradient(VT const &x) {
    return (ub - x).cwiseInverse() - (x - lb).cwiseInverse();
  }

  VT hessian(VT const &x) {
    VT d = ((ub - x).cwiseProduct((ub - x))).cwiseInverse() +
           ((x - lb).cwiseProduct((x - lb))).cwiseInverse();
    return d + extraHessian;
  }
  VT tensor(VT const &x) {
    VT d = 2 * (((ub - x).cwiseProduct((ub - x))).cwiseProduct((ub - x)))
                   .cwiseInverse() -
           2 * (((x - lb).cwiseProduct((x - lb))).cwiseProduct((x - lb)))
                   .cwiseInverse();
    return d;
  }
  VT quadratic_form_gradient(VT const &x, VT const &u) {
    // Output the -grad of u' (hess phi(x)) u.

    return (u.cwiseProduct(u)).cwiseProduct(tensor(x));
  }
  NT step_size(VT const &x, VT const &v) {
    // Output the maximum step size from x with direction v.

    // check positive direction
    VT temp = (v.array() > 0).select((ub - x).cwiseQuotient(v), max_step);
    NT t1 = temp.minCoeff();

    // check negative direction
    temp = (v.array() < 0).select((lb - x).cwiseQuotient(v), max_step);
    NT t2 = temp.minCoeff();

    return std::min(t1, t2);
  }
  VT boundary_distance(VT const &x) {
    // Output the distance of x with its closest boundary for each
    // coordinate

    return ((x - lb).cwiseMin(ub - x)).cwiseAbs();
  }

  bool feasible(VT const &x) {
    return (x.array() > lb.array() && x.array() < ub.array()).all();
  }

  std::pair<VT, VT> analytic_center_oracle(VT const &x) {
    //[~, g, h] = o.f_oracle(x);
    VT g = VT::Zero(n, 1);
    VT h = VT::Zero(n, 1);
    return std::make_pair(g + gradient(x), h + hessian(x));
  }

  std::pair<VT, VT> lewis_center_oracle(VT const &x, VT const &w) {
    //   [~, g, h] = o.f_oracle(x);
    VT g = VT::Zero(n, 1);
    VT h = VT::Zero(n, 1);
    return std::make_pair(g + w.cwiseProduct(gradient(x)),
                          h + w.cwiseProduct(hessian(x)));
  }

  std::pair<VT, VT> boundary(VT const &x) {
    // Output the normal at the boundary around x for each barrier.
    // Assume: only 1 vector is given

    VT A = VT::Ones(x.rows(), 1);

    VT b = ub;

    b = (x.array() < center.array()).select(-lb, b);

    A = (x.array() < center.array()).select(-A, A);

    return std::make_pair(A, b);
  }
};
#endif
