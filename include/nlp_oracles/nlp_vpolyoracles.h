/*
GeomScale Project

Copyright (c) 2020
  Vissarion Fisikopoulos
  Apostolos Chalkis
  Elias Tsigaridas
  Marios Papachristou

Contributed and/or modified by Marios Papachristou,
as part of Google Summer of Code 2020 program.

VolEsti is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

VolEsti is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

See the file COPYING.LESSER for the text of the GNU Lesser General
Public License.  If you did not receive this file along with HeaDDaCHe,
see <http://www.gnu.org/licenses/>.

  This file contains the implementation of boundary oracles between an
  d-dimensional V-polytope given by V with a curve p(t) = sum c_j phi_j(t)
  where c_j are d-dimensional coefficients and phi_j(t) are basis functions
  (e.g. polynomials, rational functions, splines). The boundary oracle
  returns one root (in general multiple roots exist) assuming that for t >= t_p
  the curve p(t) penetrates the V-polytope. The problem reduces to the following
  non-linear optimization problem

  max t

  subject to
              t >= 0
              lambda_i >= 0                                                   i \in [V]
              \sum_{i = 1}^m lambda_i = 1
              \sum_{i = 1}^m lambda_i v_i - \sum_{i = 1}^M c_j phi_j(t) = 0

  We use interior-point methods to solve the non-linear optimization program
  using COIN-OR ipopt and the ifopt interface for Eigen + ipopt.

*/

#ifndef NLP_VPOLYORACLES_H
#define NLP_VPOLYORACLES_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>

#define TOL 1e-7

using namespace ifopt;

// Define the variable t we use in the optimization
template <typename VT, typename NT>
class VPolyOracleVariableT : public VariableSet {
public:
  NT t, tb;

  VPolyOracleVariableT(NT t_prev, NT tb_=NT(0)): VariableSet(1, "t"), t(t_prev), tb(tb_) {};

  void SetVariables(const VT& T) override {
    t = T(0);
  }

  VectorXd GetValues() const override {
    VectorXd T(1);
    T(0) = t;
    return T;
  }

  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    bounds.at(0) = Bounds(tb, (NT) inf);
    return bounds;
  };

};

// Define the variable t we use in the optimization
template <typename VT, typename NT>
class VPolyOracleVariableLambdas : public VariableSet {
public:
  VectorXd lambdas;
  int m;

  VPolyOracleVariableLambdas(int m_): m(m_), VariableSet(m, "lambdas") {
    lambdas = VectorXd(m_);
    for (int i = 0; i < m_; i++) lambdas(i) = NT(0);
  };

  void SetVariables(const VectorXd& lambdas_) override {
    lambdas = lambdas_;
  }

  VectorXd GetValues() const override {
    return lambdas;
  }

  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    for (int i = 0; i < m; i++) bounds.at(i) = Bounds(NT(0), NT(inf));
    return bounds;
  };

};

// Define the cost function f(t) = t (ipopt takes minimization so it is -t)
template <typename VT, typename NT>
class VPolyOracleCost : public CostTerm {
public:
  int m;

  VPolyOracleCost(int m_) : CostTerm("h_poly_cost"), m(m_) {};

  NT GetCost() const override {
    VectorXd T = GetVariables()->GetComponent("t")->GetValues();
    return - T(0);
  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac) const override {
    if (var_set == "t") jac.coeffRef(0, 0) = (NT) (-1.0);
    if (var_set == "lambdas") {
      for (int i = 0; i < m; i++) jac.coeffRef(0, i+1) = NT(0);
    }
  }

};

template <typename VT, typename NT>
class VPolyoracleFeasibilityLambdas : public ConstraintSet {
public:
  int m;

  VPolyoracleFeasibilityLambdas(int m_) : ConstraintSet(1, "lambdas_simplex"), m(m_) {};

  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    bounds.at(0) = Bounds(NT(1.0), NT(1.0));
    return bounds;
  }

  VectorXd GetValues() const override {
    VectorXd lambdas = GetVariables()->GetComponent("lambdas")->GetValues();
    VectorXd sum(1);
    sum(0) = NT(0);
    for (int i = 0; i < m; i++) sum =+ lambdas(i);
    return sum;
  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override {
    if (var_set == "t") jac_block.coeffRef(0, 0) = (NT) (0);
    if (var_set == "lambdas") {
      for (int i = 0; i < m; i++) jac_block.coeffRef(0, i+1) = NT(1.0);
    }
  }

};

template <typename MT, typename VT, typename NT, typename Point, class bfunc>
class VPolyOracleFeasibilityCurve : public ConstraintSet {
public:
  int m; // number of lambdas
  int M; // number of coefficients
  int d_; // dimension
  std::vector<Point> &coeffs;
  MT &V;
  NT t0;

  bfunc phi, grad_phi;

  VPolyOracleFeasibilityCurve(MT &V_, std::vector<Point> &coeffs_, NT t0_, bfunc basis, bfunc basis_grad) :
    V(V_), coeffs(coeffs_), t0(t0_), phi(basis), grad_phi(basis_grad), ConstraintSet(V.cols() ,"curve_feasibility") {
      m = V.rows();
      M = (int) (coeffs.size());
      d_ = V.cols();
    }

  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    for (int i = 0; i < d_; i++) bounds.at(i) = Bounds(NT(0), NT(0));
    return bounds;
  }

  VectorXd GetValues() const override {
    VectorXd values(d_);
    VectorXd lambdas = GetVariables()->GetComponent("lambdas")->GetValues();
    NT t = GetVariables()->GetComponent("t")->GetValues()(0);

    for (int i = 0; i < d_; i++) {
      values(i) = NT(0);

      for (int j = 0; j < m; j++) {
        values(i) += lambdas(j) * V(j, i);
      }

      for (int j = 0; j < M; j++) {
        values(i) -= phi(t, t0, j, M) * coeffs[j][i];
      }

    }

    return values;

  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override {
    if (var_set == "t") {
      NT t = GetVariables()->GetComponent("t")->GetValues()(0);
      for (int i = 0; i < d_; i++) {
          // Derivative wrt to t
          jac_block.coeffRef(i, 0) = NT(0);

          for (int j = 0; j < M; j++) {
            jac_block.coeffRef(i, 0) -= grad_phi(t, t0, j, M) * coeffs[j][i];
          }
      }
    }

    if (var_set == "lambdas") {
      for (int i = 0; i < d_; i++) {

        for (int j = 0; j < m; j++) {
            jac_block.coeffRef(i, j+1) = V(j, i);
        }
      }
    }
  }

};

template <typename MT, typename VT, typename Point, typename NT, class bfunc>
std::pair<NT, Point> curve_intersect_vpoly_ipopt_helper(NT t_prev, NT t0, MT &V, std::vector<Point> &coeffs, bfunc phi, bfunc grad_phi) {

  Problem nlp;

  int m = V.rows();

  std::shared_ptr<VPolyOracleVariableT<VT, NT>> vpolyoraclevariable_t (new VPolyOracleVariableT<VT, NT>(t_prev, t0));
  std::shared_ptr<VPolyOracleVariableLambdas<VT, NT>> vpolyoraclevariable_lambdas (new VPolyOracleVariableLambdas<VT, NT>(m));

  nlp.AddVariableSet(vpolyoraclevariable_t);
  nlp.AddVariableSet(vpolyoraclevariable_lambdas);

  std::shared_ptr<VPolyOracleCost<VT, NT>> vpolyoraclecost (new VPolyOracleCost<VT, NT>());

  nlp.AddCostSet(vpolyoraclecost);

  std::shared_ptr<VPolyoracleFeasibilityLambdas<VT, NT>> vpolyoraclefeasibility_lambdas (new VPolyoracleFeasibilityLambdas<VT, NT>(m));
  std::shared_ptr<VPolyOracleFeasibilityCurve<MT, VT, NT, Point, bfunc>> vpolyoraclefeasibility_curve (new VPolyOracleFeasibilityCurve<MT, VT, NT, Point, bfunc>(V, coeffs, t0, phi, grad_phi));

  nlp.AddConstraintSet(vpolyoraclevariable_lambdas);
  nlp.AddConstraintSet(vpolyoraclefeasibility_curve);

  IpoptSolver ipopt;
  ipopt.SetOption("linear_solver", "mumps");
  ipopt.SetOption("jacobian_approximation", "exact");
  ipopt.SetOption("tol", TOL);
  ipopt.SetOption("print_level", 0);
  ipopt.SetOption("sb", "yes");

  ipopt.Solve(nlp);

  NT t = nlp.GetOptVariables()->GetValues()(0);

  Point p(coeffs[0].dimension());

  for (unsigned int i = 0; i < coeffs.size(); i++) {
    p += phi(t, t0, i, coeffs.size()) * coeffs[i];
  }

  return std::make_pair(t, p);

}


#endif
