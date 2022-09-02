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
#include "PackedCSparse/PackedChol.h"
#include <utility>

template <typename Polytope, typename Point, int simdLen>
class Hamiltonian
{
  using VT = typename Polytope::VT;
  using IVT =Eigen::Array<double,Eigen::Dynamic,1>;
  using BVT = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
  using NT = typename Polytope::NT;
  using MT = typename Polytope::MT;
  using Tx = FloatArray<double, simdLen>;
  using CholObj = PackedChol<simdLen, int>;
  using Opts = typename Polytope::Opts;
  using pts = std::vector<MT>;
  using Barrier = TwoSidedBarrier<Point>;
  using WeightedBarrier = WeightedTwoSidedBarrier<Point>;

public:
  bool prepared = false;
  bool forceUpdate = true;  //Update function oracle temporary varibles
  Polytope &P;
  MT hess;
  bool dUDx_empty = true;
  MT last_dUdx;
  CholObj solver;
  pts xs;
  MT x;
  MT dfx;
  MT lsc;
  VT fx;
  int n;
  int m;
  int num_runs=0;
  Barrier *barrier;
  WeightedBarrier *weighted_barrier;
  Opts &options;
  Hamiltonian(Polytope &boundaries)
      : P(boundaries), solver(CholObj(boundaries.Asp)),
        options(boundaries.options)
  {
    n = P.dimension();
    m = P.equations();
    x = MT::Zero(n,simdLen);
    xs = {x,x};
    lsc = MT::Zero(n,simdLen);
    solver.accuracyThreshold=options.solver_accuracy_threashold;
    if (options.DynamicWeight)
    {
      weighted_barrier =
          new WeightedBarrier(P.barrier.lb, P.barrier.ub, P.w_center);
      weighted_barrier->extraHessian,resize(n,simdLen);
      weighted_barrier->extraHessian=MT::Ones(n,simdLen)*options.regularization_factor;
    }
    barrier = &P.barrier;
    barrier->extraHessian,resize(n,simdLen);
    barrier->extraHessian=MT::Ones(n,simdLen)*options.regularization_factor;
  }

  // Compute H(x,v)
  VT hamiltonian(MT x, MT v)
  {
    prepare({x, v});
    pts pd = DK({x, v});
    VT K = 0.5 * (v.cwiseProduct(pd[0])).colwise().sum();
    NT logdet=(NT)solver.logdet();
    VT U = ((hess.array()).log()).colwise().sum();
    U = (U +logdet*VT::Ones(simdLen))*0.5 + fx;
    VT E = U + K;
    return E;
  }
  // Helper is nan function for vectors
  template <typename MatrixType>
  IVT is_not_nan(MatrixType x)
  {
    IVT result=IVT::Ones(1,x.cols());
    for (int i = 0; i < x.rows(); i++)
    {
      for (int j = 0; j < x.cols(); j++)
      {
        if(std::isnan(x(i,j))){
          result(j)=0;
        }
      }
    }
    return result;
  }
  // Test if the values of x and v are valid and if x is feasible
  VT feasible(MT x, MT v)
  {
    VT feasible_coordinate = VT::Ones(x.cols());
    if (options.DynamicWeight)
    {
      feasible_coordinate = weighted_barrier->feasible(x);
    }
    else
    {
      feasible_coordinate = barrier->feasible(x);
    }
    VT r= feasible_coordinate.cwiseProduct((is_not_nan(x) *is_not_nan(v)).matrix());
    return r;
  }
  // prepare the solver weighted by the hessian
  void prepare(pts const &xs)
  {
    move(xs);
    if (!prepared)
    {
      MT Hinv = hess.cwiseInverse();
      solver.decompose((Tx *)Hinv.data());
    }
    dUDx_empty = true;
    prepared = true;
  }
  // Computation of the partial derivatives of the K term
  pts DK(pts const &x_bar)
  {
    MT x = x_bar[0];
    MT v = x_bar[1];
    move(x_bar);
    MT invHessV = v.cwiseQuotient(hess);
    MT input_vector = P.Asp * invHessV;
    MT out_vector = MT::Zero(m,simdLen);
    solver.solve((Tx *)input_vector.data(), (Tx *)out_vector.data());
    MT dKdv =
        invHessV - (P.Asp.transpose() * out_vector).cwiseQuotient(hess);

    MT dKdx = MT::Zero(n,simdLen);
    if (options.DynamicWeight)
    {
      dKdx =
          weighted_barrier->quadratic_form_gradient(x, dKdv) /
          2;
    }
    else
    {
      dKdx = barrier->quadratic_form_gradient(x, dKdv) /
                   2;
    }

    return {dKdv, dKdx};
  }
  // Approximate computation of the partial derivatives of the K term
  pts approxDK(pts const &x_bar, VT &nu)
  {
    MT x = x_bar[0];
    MT v = x_bar[1];
    move(x_bar);
    MT dUdv_b = P.Asp * (v - P.Asp.transpose() * nu).cwiseQuotient(hess);
    MT out_solver = MT(nu.rows(), nu.cols());
    solver.solve((Tx *)dUdv_b.data(), (Tx *)out_solver.data());
    nu = nu + out_solver;
    MT dKdv = (v - P.Asp.transpose() * nu).cwiseQuotient(hess);
    MT dKdx = MT::Zero(n,simdLen);
    if (options.DynamicWeight)
    {
      dKdx =
          weighted_barrier->quadratic_form_gradient(x, dKdv) /
          2;
    }
    else
    {
      dKdx = barrier->quadratic_form_gradient(x, dKdv) /
                   2;
    }
    return {dKdv, dKdx};
  }
  inline MT operator+(const MT &M, const VT v) {
    return M.colwise()+v;
  }
  inline MT operator-(const MT &M, const VT v) {
    return M.colwise()-v;
  }
  inline MT operator-( const VT v,const MT &M) {
    return (-M).colwise()+v;
  }
  // Compute the partial derivatives of one term
  // This is only dependent on x and so DU/Dv=0
  pts DU(pts const &x_bar)
  {
    MT x = x_bar[0];
    move(x_bar);
    if (!prepared || dUDx_empty)
    {
      prepare(x_bar);
      solver.leverageScoreComplement((Tx *)lsc.data());
      if(num_runs<10){
        std::cerr<<"--------lsc----------\n";
        std::cerr<<lsc<<"\n";
        std::cerr<<"------x--------------\n";
        std::cerr<<x<<"\n";
        std::cerr<<"------b-x.col--------------\n";
        std::cerr<<barrier->ub-x.colwise()<<"\n";
        std::cerr<<"------b-x.col^2--------------\n";
        std::cerr<<(barrier->ub-x.colwise()).cwiseProduct(barrier->ub-x.colwise())<<"\n";
        std::cerr<<"------b-x.col^2--------------\n";
        std::cerr<<(barrier->ub-x).cwiseProduct(barrier->ub-x).cwiseInverse()+(x-barrier->lb).cwiseProduct(x-barrier->lb).cwiseInverse()+<<"\n";
        std::cerr<<"------tensor--------------\n";
        std::cerr<< barrier->tensor(x)<<"\n";
        std::cerr<<"---------------end------------\n";
      }
      if (options.DynamicWeight)
      {
        last_dUdx = (weighted_barrier->tensor(x).cwiseProduct(lsc))
                               .cwiseQuotient(2 * hess) +
                          dfx;
      }
      else
      {
        last_dUdx =
            (barrier->tensor(x).cwiseProduct(lsc)).cwiseQuotient(2 * hess) +
            dfx;
      }
      dUDx_empty = false;
    }

    return {MT::Zero(n,simdLen), -last_dUdx};
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
    x = xs[0];
    MT h;
    std::tie(fx, dfx, h) = P.f_oracle(x);
    num_runs++;
    if(num_runs<10){
      std::cerr<<"--------fx----------\n";
      std::cerr<<fx<<"\n";
      std::cerr<<"--------dfx----------\n";
      std::cerr<<dfx<<"\n";
      std::cerr<<"--------h----------\n";
      std::cerr<<h<<"\n";
      std::cerr<<"--------end----------\n";

    }
    if (options.DynamicWeight)
    {
      hess = weighted_barrier->hessian(x) + h;
    }
    else
    {
      hess = barrier->hessian(x) + h;
    }
    if(num_runs<10){
      std::cerr<<"----------barrier->hess-------------\n";
      std::cerr << barrier->hessian(x) << '\n';
      std::cerr<<"--------hess----------\n";
      std::cerr<<hess<<"\n";
      std::cerr<<"--------end----------\n";
    }
    forceUpdate = false;
    prepared = false;
  }
  // Project x to the polytope
  void project(pts &xs)
  {
    move(xs);
    MT x = xs[0];
    int m = P.Asp.rows();
    MT out_vector = MT(m,simdLen);
    MT in_vector = P.b - P.Asp * x;
    solver.solve((Tx *)in_vector.data(), (Tx *)out_vector.data());
    out_vector = P.Asp.transpose() * out_vector;
    xs[0] = xs[0] + (out_vector).cwiseQuotient(hess);
  }
  // Get the inner product of x and ds weighted by the hessian
  VT x_norm(pts const &xs, pts const &dx)
  {
    move(xs);
    MT dx_x = dx[0];
    MT r = (dx_x.cwiseProduct(dx_x)).cwiseProduct(hess);
    return r.colwise().sum();
  }
};
#endif
