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
*/

#ifndef NLP_HPOLYORACLES_H
#define NLP_HPOLYORACLES_H

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

template <typename VT, typename NT>
class HPolyOracleVariables : public VariableSet {
public:
  NT t, tb;

  HPolyOracleVariables(NT t_prev): VariableSet(1, "t"), t(t_prev), tb(NT(0)) {};

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

template <typename VT, typename NT>
class HPolyOracleCost : public CostTerm {
public:
  HPolyOracleCost() : CostTerm("h_poly_cost") {};

  NT GetCost() const override {
    VectorXd T = GetVariables()->GetComponent("t")->GetValues();
    return - T(0);
  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac) const override {
    if (var_set == "t") jac.coeffRef(0, 0) = (NT) (-1.0);
  }

};

template <typename MT, typename VT, typename NT, class bfunc>
class HPolyOracleFeasibility : public ConstraintSet {
public:
  MT &C;
  VT &b;
  bfunc phi, grad_phi;
  NT t0;
  int m, M;

  HPolyOracleFeasibility(MT &C_, VT &b_, NT t0_, bfunc basis, bfunc basis_grad) :
    C(C_), b(b_), t0(t0_), phi(basis), grad_phi(basis_grad), ConstraintSet(C_.rows(), "h_poly_feasibility") {
      m = C_.rows();
      M = C_.cols();
    };

  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    for (int i = 0; i < m; i++) {
      bounds.at(i) = Bounds((NT) (-inf), b(i));
    }
    return bounds;
  }

  VectorXd GetValues() const override {
    VectorXd T = GetVariables()->GetComponent("t")->GetValues();
    NT t = T(0);
    VectorXd phis;
    phis.resize(M);

    for (int i = 0; i < M; i++) {
      phis(i) = (NT) (phi(t, t0, i, M));
    }
    return C * phis;
  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override {

    if (var_set == "t") {
      VectorXd T = GetVariables()->GetComponent("t")->GetValues();
      NT t = T(0);
      NT temp;

      for (int i = 0; i < m; i++) {
        jac_block.coeffRef(i, 0) = NT(0);
        for (int j = 0; j < M; j++) {
          temp = grad_phi(t, t0, j, M);
          if (std::isinf(temp)) continue;
          else jac_block.coeffRef(i, 0) += temp;
        }
      }
    }
  }


};

template <typename MT, typename VT, typename Point, typename NT, class bfunc>
std::tuple<NT, Point, int> curve_intersect_ipopt_helper(NT t_prev, NT t0, MT &A, VT &b, std::vector<Point> &coeffs, bfunc phi, bfunc grad_phi)
{

  Problem nlp;

  MT C, C_tmp;
  C_tmp.resize(coeffs[0].dimension(), coeffs.size());


  for (int i = 0; i < coeffs.size(); i++) {
    C_tmp.col(i) = coeffs[i].getCoefficients();
  }

  // C_tmp: dimension x num_coeffs
  // A: constraints x dimension
  // C: constraints x num_coeffs
  C = A * C_tmp;

  std::shared_ptr<HPolyOracleVariables<VT, NT>> hpolyoraclevariables (new HPolyOracleVariables<VT, NT>(t_prev));
  std::shared_ptr<HPolyOracleCost<VT, NT>> hpolyoraclecost (new HPolyOracleCost<VT, NT>());
  std::shared_ptr<HPolyOracleFeasibility<MT, VT, NT, bfunc>> hpolyoraclefeasibility (new HPolyOracleFeasibility<MT, VT, NT, bfunc>(C, b, t0, phi, grad_phi));

  nlp.AddVariableSet  (hpolyoraclevariables);
  nlp.AddCostSet      (hpolyoraclecost);
  nlp.AddConstraintSet(hpolyoraclefeasibility);
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

  const NT* b_data = b.data();


  for (int i = 0; i < A.rows(); i++) {
    if (*b_data == 0 && std::abs(*b_data - A.row(i) * p.getCoefficients()) < NT(1000 * TOL)) {
      return std::make_tuple(t, p, i);
    }
    else if (*b_data != 0 && std::abs(*b_data - A.row(i) * p.getCoefficients()) / *b_data < NT(1000 * TOL)) {
      return std::make_tuple(t, p, i);
    }
    else b_data++;
  }

  return std::make_tuple(t, p, -1);
}


#endif
