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
#ifndef ANALYTIC_CENTER_H
#define ANALYTIC_CENTER_H
#include "Eigen/Eigen"
#include "PackedCSparse/PackedChol.h"
#include "preprocess/crhmc/crhmc_utils.h"
#include "preprocess/crhmc/opts.h"
#include "sos/barriers/TwoSidedBarrier.h"
#include <fstream>
#include <iostream>

#include <vector>
#ifndef SIMD_LEN
#define SIMD_LEN 0
#endif
const size_t chol_k2 = (SIMD_LEN == 0) ? 1 : SIMD_LEN;

using NT=double;
using MT=Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
using VT=Eigen::Matrix<NT, Eigen::Dynamic, 1>;
using SpMat= Eigen::SparseMatrix<NT>;
using CholObj=PackedChol<chol_k2, int>;
using Triple=Eigen::Triplet<double>;
using Barrier=TwoSidedBarrier<NT>;
using Tx = FloatArray<double, chol_k2>;
using Opts=opts<NT>;

std::tuple<VT, SpMat, VT> analytic_center(SpMat const &A, VT const &b,
                                          Barrier &f, Opts const &options,
                                          VT x = VT::Zero(0, 1)) {
  // initial conditions
  int n = A.cols();
  int m = A.rows();
  if (x.rows() == 0 || !f.feasible(x)) {
    x = f.center;
  }

  VT lambda = VT::Zero(n, 1);
  int fullStep = 0;
  NT tConst = 0;
  NT primalErr = std::numeric_limits<NT>::infinity();
  NT dualErr = std::numeric_limits<NT>::infinity();
  NT primalErrMin = std::numeric_limits<NT>::infinity();
  NT primalFactor = 1;
  NT dualFactor = 1 + b.norm();
  std::vector<int> idx;

  CholObj solver = CholObj(A);

  for (int iter = 0; iter < options.ipmMaxIter; iter++) {
    std::pair<VT, VT> pair_analytic_oracle = f.analytic_center_oracle(x);
    VT grad = pair_analytic_oracle.first;
    VT hess = pair_analytic_oracle.second;

    // compute the residual
    VT rx = lambda - grad;
    VT rs = b - A * x;

    // check stagnation
    NT primalErrMin = std::min(primalErr, primalErrMin);
    NT primalErr = rx.norm() / primalFactor;
    NT dualErrLast = dualErr;
    NT dualErr = rs.norm() / dualFactor;
    bool feasible = f.feasible(x);
    if ((dualErr > (1 - 0.9 * tConst) * dualErrLast) ||
        (primalErr > 10 * primalErrMin) || !feasible) {
      VT dist = f.boundary_distance(x);
      NT th = options.ipmDistanceTol;
      visit_lambda(dist, [&idx, th](double v, int i, int j) {
        if (v < th)
          idx.push_back(i);
      });

      // idx = find(dist < ipmDistanceTol);
      if (idx.size() > 0) {
        break;
      }
    }

    // compute the step direction
    VT Hinv = hess.cwiseInverse();
    solver.decompose(Hinv.data());
    VT out(m, 1);
    solver.solve(rs.data(), out.data());
    VT dr1 = A.transpose() * out;
    VT in = A * Hinv.cwiseProduct(rx);
    solver.solve(in.data(), out.data());

    VT dr2 = A.transpose() * out;
    VT dx1 = Hinv.cwiseProduct(dr1);
    VT dx2 = Hinv.cwiseProduct(rx - dr2);

    // compute the step size
    VT dx = dx1 + dx2;
    NT tGrad = std::min(f.step_size(x, dx), 1.0);
    dx = dx1 + tGrad * dx2;
    NT tConst = std::min(0.99 * f.step_size(x, dx), 1.0);
    tGrad = tGrad * tConst;

    // make the step
    x = x + tConst * dx;
    lambda = lambda - dr2;

    if (!f.feasible(x)) {
      break;
    }

    if (tGrad == 1) {
      fullStep = fullStep + 1;
      if (fullStep > log(dualErr / options.ipmDualTol) && fullStep > 8) {
        break;
      }
    } else {
      fullStep = 0;
    }
  }
  SpMat C;
  VT d;
  if (idx.size() == 0) {
    VT dist = f.boundary_distance(x);
    NT th = options.ipmDistanceTol;
    visit_lambda(dist, [&idx, th](double v, int i, int j) {
      if (v < th)
        idx.push_back(i);
    });
    // idx = find(dist < ipmDistanceTol);
  }

  if (idx.size() > 0) {
    std::pair<VT, VT> pboundary = f.boundary(x);
    VT A_ = pboundary.first;
    VT b_ = pboundary.second;
    A_ = A_(idx);
    std::vector<Triple> sparseIdx;
    for (int i = 0; i < idx.size(); i++) {
      sparseIdx.push_back(Triple(i, i, A_(i)));
    }
    C.setFromTriplets(sparseIdx.begin(), sparseIdx.end());
    d = b_(idx);
  } else {
    C = MT::Zero(0, n).sparseView();
    d = VT::Zero(0, 1);
  }
  return std::make_tuple(x, C, d);
}
#endif
