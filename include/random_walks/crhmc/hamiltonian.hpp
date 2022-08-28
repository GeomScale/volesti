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
#include "sos/barriers/TwoSidedBarrier.h"
#include "sos/barriers/WeightedTwoSidedBarrier.h"
#include <utility>

template <typename Polytope, typename Point>
class Hamiltonian
{
  using VT = typename Polytope::VT;
  using NT = typename Polytope::NT;
  using MT = typename Polytope::MT;
  using Tx = typename Polytope::Tx;
  using CholObj = typename Polytope::CholObj;
  using Opts = typename Polytope::Opts;

  using pts = std::vector<Point>;
  using Barrier = TwoSidedBarrier<Point>;
  using WeightedBarrier = WeightedTwoSidedBarrier<Point>;

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
  VT lsc;
  NT fx = 0;
  int n;
  int m;
  Barrier *barrier;
  WeightedBarrier *weighted_barrier;
  Opts &options;
  Hamiltonian(Polytope &boundaries)
      : P(boundaries), solver(CholObj(boundaries.Asp)),
        options(boundaries.options)
  {
    n = P.dimension();
    m = P.equations();
    x = VT::Zero(n);
    xs = {Point(n), Point(n)};
    lsc = VT::Zero(n);
    solver.accuracyThreshold=options.solver_accuracy_threashold;
    if (options.DynamicWeight)
    {
      weighted_barrier =
          new WeightedBarrier(P.barrier.lb, P.barrier.ub, P.w_center);
    }
    barrier = &P.barrier;
  }

  // Compute H(x,v)
  NT hamiltonian(Point x, Point v)
  {
    prepare({x, v});
    pts pd = DK({x, v});
    NT K = 0.5 * v.dot(pd[0]);
    NT U = 0.5 * (solver.logdet() + ((hess.array()).log()).sum());
    U = U + fx;
    NT E = U + K;
    return E;
  }
  // Helper is nan function for vectors
  template <typename MatrixType>
  bool isnan(MatrixType x)
  {
    for (int i = 0; i < x.rows(); i++)
    {
      for (int j = 0; j < x.cols(); j++)
      {
        if (std::isnan(x(i, j)))
        {
          return true;
        }
      }
    }
    return false;
  }
  // Test if the values of x and v are valid and if x is feasible
  NT feasible(VT x, VT v)
  {
    bool feasible_coordinate = true;
    if (options.DynamicWeight)
    {
      feasible_coordinate = weighted_barrier->feasible(x);
    }
    else
    {
      feasible_coordinate = barrier->feasible(x);
    }
    bool r = !isnan(x) && !isnan(v) && feasible_coordinate;
    if (r)
    {
      return 1;
    }
    return 0;
  }
  // prepare the solver weighted by the hessian
  void prepare(pts const &xs)
  {
    move(xs);
    if (!prepared)
    {
      VT Hinv = hess.cwiseInverse();
      solver.decompose((Tx *)Hinv.data());
    }
    dUDx_empty = true;
    prepared = true;
  }
  // Computation of the partial derivatives of the K term
  pts DK(pts const &x_bar)
  {
    VT x = x_bar[0].getCoefficients();
    VT v = x_bar[1].getCoefficients();
    move(x_bar);
    VT invHessV = v.cwiseQuotient(hess);
    VT input_vector = P.Asp * invHessV;
    VT out_vector = VT::Zero(m);
    solver.solve((Tx *)input_vector.data(), (Tx *)out_vector.data());
    Point dKdv =
        Point(invHessV - (P.Asp.transpose() * out_vector).cwiseQuotient(hess));

    Point dKdx = Point(n);
    if (options.DynamicWeight)
    {
      dKdx = Point(
          weighted_barrier->quadratic_form_gradient(x, dKdv.getCoefficients()) /
          2);
    }
    else
    {
      dKdx = Point(barrier->quadratic_form_gradient(x, dKdv.getCoefficients()) /
                   2);
    }

    return {dKdv, dKdx};
  }
  // Approximate computation of the partial derivatives of the K term
  pts approxDK(pts const &x_bar, VT &nu)
  {
    VT x = x_bar[0].getCoefficients();
    VT v = x_bar[1].getCoefficients();
    move(x_bar);
    VT dUdv_b = P.Asp * (v - P.Asp.transpose() * nu).cwiseQuotient(hess);
    VT out_solver = VT(nu.rows(), nu.cols());
    solver.solve((Tx *)dUdv_b.data(), (Tx *)out_solver.data());
    nu = nu + out_solver;
    Point dKdv = Point((v - P.Asp.transpose() * nu).cwiseQuotient(hess));
    Point dKdx = Point(n);
    if (options.DynamicWeight)
    {
      dKdx = Point(
          weighted_barrier->quadratic_form_gradient(x, dKdv.getCoefficients()) /
          2);
    }
    else
    {
      dKdx = Point(barrier->quadratic_form_gradient(x, dKdv.getCoefficients()) /
                   2);
    }
    return {dKdv, dKdx};
  }
  // Compute the partial derivatives of one term
  // This is only dependent on x and so DU/Dv=0
  pts DU(pts const &x_bar)
  {
    VT x = x_bar[0].getCoefficients();
    move(x_bar);
    if (!prepared || dUDx_empty)
    {
      prepare(x_bar);
      solver.leverageScoreComplement((Tx *)lsc.data());

      if (options.DynamicWeight)
      {
        last_dUdx = Point(-(weighted_barrier->tensor(x).cwiseProduct(lsc))
                               .cwiseQuotient(2 * hess) -
                          dfx);
      }
      else
      {
        last_dUdx = Point(
            -(barrier->tensor(x).cwiseProduct(lsc)).cwiseQuotient(2 * hess) -
            dfx);
      }
      dUDx_empty = false;
    }
    return {Point(n), last_dUdx};
  }
  // Compute the computations involving only x iff x has been changed
  // Else they are stored
  void move(pts const &y)
  {
    if (y[0] == xs[0] && !forceUpdate)
    {
      return;
    }
    xs = y;
    x = xs[0].getCoefficients();
    VT h;
    std::tie(fx, dfx, h) = P.f_oracle(x);
    if (options.DynamicWeight)
    {
      hess = weighted_barrier->hessian(x) + h;
    }
    else
    {
      hess = barrier->hessian(x) + h;
    }

    forceUpdate = false;
    prepared = false;
  }
  // Project x to the polytope
  void project(pts &xs)
  {
    move(xs);
    VT x = xs[0].getCoefficients();
    int m = P.Asp.rows();
    VT out_vector = VT(m);
    VT in_vector = P.b - P.Asp * x;
    solver.solve((Tx *)in_vector.data(), (Tx *)out_vector.data());
    out_vector = P.Asp.transpose() * out_vector;
    xs[0] = xs[0] + Point((out_vector).cwiseQuotient(hess));
  }
  // Get the inner product of x and ds weighted by the hessian
  NT x_norm(pts const &xs, pts const &dx)
  {
    move(xs);
    VT dx_x = dx[0].getCoefficients();
    VT r = (dx_x.cwiseProduct(dx_x)).cwiseProduct(hess);
    return r.sum();
  }
};
#endif
