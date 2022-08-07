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
#ifndef HAMILTONIAN_UTILS_HPP
#define HAMILTONIAN_UTILS_HPP
#include <utility>
template <typename Polytope, typename func> class Hamiltonian {
  using VT = typename Polytope::VT;
  using NT = typename Polytope::NT;
  using MT = typename Polytope::MT;
  using Tx = typename Polytope::Tx;

  typedef typename Polytope::CholObj CholObj;
  typedef std::vector<VT> pts;

public:
  bool prepared = false;
  bool forceUpdate = true;
  Polytope *P;
  VT hess;
  bool dUDx_empty = true;
  VT last_dUdx;
  CholObj solver;
  VT x;
  VT dfx;

  func F;
  int n;
  int m;

  Hamiltonian(Polytope *boundaries, func oracle)
      : P(boundaries), F(oracle), solver(CholObj(boundaries->Asp)) {
    n = P->dimension();
    m = P->equations();
    x = VT::Zero(n);
  }
  void prepare(VT const &x) {
    move(x);
    if (!prepared) {
      VT Hinv = hess.cwiseInverse();
      solver.decompose((Tx *)Hinv.data());
    }
    dUDx_empty = true;
    prepared = true;
  }
  pts DK(pts const &x_bar) {
    VT x = x_bar[0];
    VT v = x_bar[1];
    move(x);
    VT invHessV = v.cwiseQuotient(hess);
    VT input_vector = P->Asp * invHessV;
    VT out_vector = VT::Zero(m);
    solver.solve((Tx *)input_vector.data(), (Tx *)out_vector.data());
    VT dKdv = invHessV - (P->A.transpose() * out_vector).cwiseQuotient(hess);
    VT dKdx = -P->barrier.quadratic_form_gradient(x, dKdv) / 2;
    return {dKdv, -dKdx};
  }
  std::pair<pts, MT> approxDK(pts const &x_bar, MT const &nu) {
    VT x = x_bar[0];
    VT v = x_bar[1];
    move(x);
    VT dUdv_b = P->Asp * (v - P->Asp.transpose() * nu).cwiseQuotient(hess);

    VT out_solver = VT(nu.rows(), nu.cols());
    solver.solve((Tx *)dUdv_b.data(), (Tx *)out_solver.data());
    nu = nu + out_solver;

    VT dKdv = (v - P->Asp.transpose() * nu).cwiseQuotient(hess);
    VT dKdx = -P->barrier.quadratic_form_gradient(x, dKdv) / 2;
    pts result = {dUdv_b, -dKdx};
    return std::make_pair(result, nu);
  }

  pts DU(pts const &x_bar) {
    VT x = x_bar[0];
    move(x);
    if (!prepared || dUDx_empty) {
      prepare(x);
      VT lsc = VT(n, 1);
      solver.leverageScoreComplement((Tx *)lsc.data());

      last_dUdx =
          (P->barrier.tensor(x).cwiseProduct(lsc)).cwiseQuotient(2 * hess) +
          dfx;
      dUDx_empty = false;
    }
    return {VT::Zero(x.rows(), 1), -last_dUdx};
  }

  void move(VT const &y) {
    if (y.isApprox(x) || forceUpdate) {
      return;
    }

    x = y;
    dfx = F(x);
    hess = P->barrier.hessian(x);
    forceUpdate = false;
    prepared = false;
  }

  VT x_norm(VT const &x, VT const &dx) {
    move(x);
    VT r = (dx.cwiseProduct(dx)).cwiseProduct(hess);
    return r;
  }
};
#endif
