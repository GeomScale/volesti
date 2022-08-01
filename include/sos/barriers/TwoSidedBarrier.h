// The log barrier for the domain {lu <= x <= ub}:
//	phi(x) = - sum log(x - lb) - sum log(ub - x).
#ifndef TWOSIDEDBARIER_H
#define TWOSIDEDBARIER_H

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include <vector>

template <typename Point> class TwoSidedBarrier {

  typedef double NT;
  typedef Cartesian<NT> Kernel;
  typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
  typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
  typedef Eigen::SparseMatrix<NT> SpMat;

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
  const NT inf = std::numeric_limits<NT>::infinity();

  void set_bound(VT const &_lb, VT const &_ub) {

    lb = _lb;
    ub = _ub;
    n = lb.rows();
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
  }

  VT gradient(VT const &x) {
    return (ub - x).cwiseInverse() - (x - lb).cwiseInverse();
  }

  VT hessian(VT const &x) {
    VT d = ((ub - x).cwiseProduct((ub - x))).cwiseInverse() +
           ((x - lb).cwiseProduct((x - lb))).cwiseInverse();
    return d;
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