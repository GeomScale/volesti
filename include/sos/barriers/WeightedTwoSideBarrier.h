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
#ifndef WEIGHTEDTWOSIDEDBARIER_H
#define WEIGHTEDTWOSIDEDBARIER_H

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include <vector>

template <typename Point> class WeightedTwoSidedBarrier:public TwoSidedBarrier<Point> {

  using NT = Point::FP;
  using Kernel = Cartesian<NT>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using SpMat = Eigen::SparseMatrix<NT>;

public:

  VT w;

  WeightedTwoSidedBarrier(VT const &_lb, VT const &_ub, int _vdim = 1) {
    set_bound(_lb, _ub);
    vdim = _vdim;
  }
  WeightedTwoSidedBarrier() { vdim = 1; }

  VT gradient(VT const &x) {
    return w.cwiseQuotient(ub - x) - w.cwiseQuotient(x - lb);
  }

  VT hessian(VT const &x) {
    VT d = w.cwiseQuotient((ub - x).cwiseProduct((ub - x))) +
           w.cwiseQuotient((x - lb).cwiseProduct((x - lb)));
    return d;
  }
  VT tensor(VT const &x) {
    VT d = 2 * w.cwiseQuotient(((ub - x).cwiseProduct((ub - x))).cwiseProduct((ub - x)))
                    -
           2 * w.cwiseQuotient(((x - lb).cwiseProduct((x - lb))).cwiseProduct((x - lb)));
    return d;
  }


};
#endif
