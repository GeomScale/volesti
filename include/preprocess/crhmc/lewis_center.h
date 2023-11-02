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
#ifndef LEWIS_CENTER_H
#define LEWIS_CENTER_H
#include "Eigen/Eigen"
#include "PackedCSparse/PackedChol.h"
#include "preprocess/crhmc/crhmc_utils.h"
#include "preprocess/crhmc/opts.h"
#include "preprocess/crhmc/two_sided_barrier.h"
#include <fstream>
#include <iostream>
#include <vector>
#ifndef SIMD_LEN
#define SIMD_LEN 0
#endif
const size_t chol_k3 = (SIMD_LEN == 0) ? 1 : SIMD_LEN;

/*This function computes the Lewis center of the polytope*/
//And detects additional constraint that need to be added
// x - It outputs the minimizer of min f(x) subjects to {Ax=b}
// w - Output weights that correspond to the Lewis center, they are gone be used in the sampler to reduce the conditon number
// C - detected constraint matrix
//     If the domain ({Ax=b} intersect dom(f)) is not full dimensional in {Ax=b}
//     because of the dom(f), the algorithm will detect the collapsed dimension
//     and output the detected constraint C x = d
// d - detected constraint vector
template <typename Polytope, typename SpMat,typename Opts, typename MT, typename VT, typename NT>
std::tuple<VT, SpMat, VT, VT> lewis_center(SpMat const &A, VT const &b, Polytope &f, Opts const &options, VT x = VT::Zero(0, 1))
{
  using CholObj = typename Polytope::CholObj;
  using Triple = typename Polytope::Triple;
  using Tx = typename Polytope::Tx;
  NT epsilon = 1e-8;
  // initial conditions
  int n = A.cols();
  int m = A.rows();
  //If it is given use starting point
  if (x.rows() == 0 || !f.barrier.feasible(x))
  {
    x = f.barrier.center;
  }
  VT lambda = VT::Zero(n, 1);
  int fullStep = 0;
  NT tConst = 0;
  NT primalErr = std::numeric_limits<NT>::max();
  NT dualErr = std::numeric_limits<NT>::max();
  NT primalErrMin = std::numeric_limits<NT>::max();
  NT primalFactor = 1;
  NT dualFactor = 1 + b.norm();
  std::vector<int> idx;

  CholObj solver = CholObj(transform_format<SpMat,NT,int>(A));
  VT w = VT::Ones(n, 1);
  VT wp = w;
  for (int iter = 0; iter < options.ipmMaxIter; iter++)
  {
    std::pair<VT, VT> pair_analytic_oracle = f.lewis_center_oracle(x, wp);
    VT grad = pair_analytic_oracle.first;
    VT hess = pair_analytic_oracle.second;

    // compute the residual
    VT rx = lambda - grad;
    VT rs = b - A * x;

    // check stagnation
    primalErrMin = std::min(primalErr, primalErrMin);
    primalErr = rx.norm() / primalFactor;
    NT dualErrLast = dualErr;
    dualErr = rs.norm() / dualFactor;
    bool feasible = f.barrier.feasible(x);
    if ((dualErr > (1 - 0.9 * tConst) * dualErrLast) ||
        (primalErr > 10 * primalErrMin) || !feasible)
    {
      VT dist = f.barrier.boundary_distance(x);
      NT th = options.ipmDistanceTol;
      visit_lambda(dist, [&idx, th](double v, int i, int j)
                   {
        if (v < th)
          idx.push_back(i); });

      if (idx.size() > 0)
      {
        break;
      }
    }

    // compute the step direction
    VT Hinv = hess.cwiseInverse();
    solver.decompose((Tx *)Hinv.data());
    VT out(m, 1);
    solver.solve((Tx *)rs.data(), (Tx *)out.data());
    VT dr1 = A.transpose() * out;
    VT in = A * Hinv.cwiseProduct(rx);
    solver.solve((Tx *)in.data(), (Tx *)out.data());

    VT dr2 = A.transpose() * out;
    VT dx1 = Hinv.cwiseProduct(dr1);
    VT dx2 = Hinv.cwiseProduct(rx - dr2);

    // compute the step size
    VT dx = dx1 + dx2;
    NT tGrad = std::min(f.barrier.step_size(x, dx), 1.0);
    dx = dx1 + tGrad * dx2;
    NT tConst = std::min(0.99 * f.barrier.step_size(x, dx), 1.0);
    tGrad = tGrad * tConst;

    // make the step
    x = x + tConst * dx;
    lambda = lambda - dr2;

    // update weight
    VT w_vector(n, 1);
    solver.leverageScoreComplement((Tx *)w_vector.data());

    VT wNew = w_vector.cwiseMax(0) + VT::Ones(n, 1) * epsilon;
    w = (w + wNew) / 2;
    wp = Eigen::pow(w.array(), 0.875).matrix();

    if (!f.barrier.feasible(x))
    {
      break;
    }

    // stop if converged
    if (tGrad == 1)
    {
      fullStep = fullStep + 1;
      if (fullStep > log(dualErr / options.ipmDualTol) &&
          fullStep > options.min_convergence_steps)
      {
        break;
      }
    }
    else
    {
      fullStep = 0;
    }
  }

  SpMat C;
  VT d;
  if (idx.size() == 0)
  {
    VT dist = f.barrier.boundary_distance(x);
    NT th = options.ipmDistanceTol;
    visit_lambda(dist, [&idx, th](double v, int i, int j)
                 {
      if (v < th)
        idx.push_back(i); });
  }

  if (idx.size() > 0)
  {
    C.resize(idx.size(), n);
    std::pair<VT, VT> pboundary = f.barrier.boundary(x);
    VT A_ = pboundary.first;
    VT b_ = pboundary.second;
    std::vector<Triple> sparseIdx;
    for (int i = 0; i < idx.size(); i++)
    {
      sparseIdx.push_back(Triple(i, idx[i], A_(idx[i])));
    }
    C.setFromTriplets(sparseIdx.begin(), sparseIdx.end());
    d.resize(idx.size(), 1);
    copy_indicies(d, b_, idx);
  }
  else
  {
    C = MT::Zero(0, n).sparseView();
    d = VT::Zero(0, 1);
  }

  return std::make_tuple(x, C, d, wp);
}
#endif
