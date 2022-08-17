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
#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP
#include <utility>
template <typename Polytope, typename Point, typename func> class Hamiltonian {
  using VT = typename Polytope::VT;
  using NT = typename Polytope::NT;
  using MT = typename Polytope::MT;
  using Tx = typename Polytope::Tx;
  using CholObj = typename Polytope::CholObj;
  using pts = std::vector<Point>;

public:
  bool prepared = false;
  bool forceUpdate = true;
  Polytope &P;
  VT hess;
  bool dUDx_empty = true;
  Point last_dUdx;
  CholObj solver;
  pts xs;
  VT x;
  VT dfx;
  NT fx = 0;
  func F;
  int n;
  int m;

  Hamiltonian(Polytope &boundaries, func oracle)
      : P(boundaries), F(oracle), solver(CholObj(boundaries.Asp)) {
    n = P.dimension();
    m = P.equations();
    x = VT::Zero(n);
    xs = {Point(n), Point(n)};
  }

  // Compute H(x,v)
  NT hamiltonian(Point x, Point v) {
    prepare({x, v});
    pts pd = DK({x, v});
    NT K = 0.5 * v.dot(pd[0]);
    NT U = 0.5 * (solver.logdet() + ((hess.array()).log()).sum());
    U = U + fx;
    NT E = U + K;
    std::cerr << "x= (" << x.getCoefficients().transpose() << ") v= ("
              << -v.getCoefficients().transpose() << ")\n";
    return E;
  }
  template <typename MatrixType> bool isnan(MatrixType x) {
    for (int i = 0; i < x.rows(); i++) {
      for (int j = 0; j < x.cols(); j++) {
        if (std::isnan(x(i, j)))
          return 1;
      }
    }
    return 0;
  }
  // Test if the values of x and v are valid and if x is feasible
  NT feasible(VT x, VT v) {

    bool r = !isnan(x) && !isnan(v) && P.barrier.feasible(x);
    if (r) {
      return 1;
    }
    return 0;
  }
  void prepare(pts const &xs) {
    move(xs);
    if (!prepared) {
      VT Hinv = hess.cwiseInverse();
      solver.decompose((Tx *)Hinv.data());
    }
    dUDx_empty = true;
    prepared = true;
  }
  pts DK(pts const &x_bar) {
    VT x = x_bar[0].getCoefficients();
    VT v = x_bar[1].getCoefficients();
    move(x_bar);
    VT invHessV = v.cwiseQuotient(hess);
    VT input_vector = P.Asp * invHessV;
    VT out_vector = VT::Zero(m);
    solver.solve((Tx *)input_vector.data(), (Tx *)out_vector.data());
    Point dKdv =
        Point(invHessV - (P.A.transpose() * out_vector).cwiseQuotient(hess));
    Point dKdx =
        Point(P.barrier.quadratic_form_gradient(x, dKdv.getCoefficients()) / 2);
    return {dKdv, dKdx};
  }
  std::pair<pts, MT> approxDK(pts const &x_bar, MT const &nu) {
    VT x = x_bar[0].getCoefficients();
    VT v = x_bar[1].getCoefficients();
    move(x_bar);
    VT dUdv_b = P.Asp * (v - P.Asp.transpose() * nu).cwiseQuotient(hess);

    VT out_solver = VT(nu.rows(), nu.cols());
    solver.solve((Tx *)dUdv_b.data(), (Tx *)out_solver.data());
    nu = nu + out_solver;

    Point dKdv = Point((v - P.Asp.transpose() * nu).cwiseQuotient(hess));
    Point dKdx =
        Point(P.barrier.quadratic_form_gradient(x, dKdv.getCoefficients()) / 2);
    pts result = {dUdv_b, dKdx};
    return std::make_pair(result, nu);
  }

  pts DU(pts const &x_bar) {
    VT x = x_bar[0].getCoefficients();
    move(x_bar);
    if (!prepared || dUDx_empty) {
      prepare(x_bar);
      VT lsc = VT(n, 1);
      solver.leverageScoreComplement((Tx *)lsc.data());

      last_dUdx = Point(
          -(P.barrier.tensor(x).cwiseProduct(lsc)).cwiseQuotient(2 * hess) -
          dfx);
      dUDx_empty = false;
    }
    return {Point(n), last_dUdx};
  }

  void move(pts const &y) {
    if (y[0] == xs[0] && !forceUpdate) {
      return;
    }
    xs = y;
    x = xs[0].getCoefficients();
    std::tie(fx, dfx, std::ignore) = P.f_oracle(x);
    hess = P.barrier.hessian(x);
    forceUpdate = false;
    prepared = false;
  }

  VT x_norm(pts const &xs, pts const &dx) {
    move(xs);
    VT dx_x = dx[0].getCoefficients();
    VT r = (dx_x.cwiseProduct(dx_x)).cwiseProduct(hess);
    return r;
  }
};
#endif
