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
  d-dimensional H-polytope P: Ax <= b with a curve p(t) = sum c_j phi_j(t)
  where c_j are d-dimensional coefficients and phi_j(t) are basis functions
  (e.g. polynomials, rational functions, splines). The boundary oracle
  returns one root (in general multiple roots exist) assuming that for t >= t_p
  the curve p(t) penetrates the H-polytope. The problem reduces to the following
  non-linear optimization problem for finding the maximum t such that p(t) penetrates
  the polytope where t lies inside [t0, t0 + eta] for some eta > 0

  max t

  subject to
              t >= t0
              t <= t0 + eta
              A p(t) <= b

  The second constraint can be rewritten as (A*C) * Phi <= b which is eventually
  the optimization problem we solve where the vector Phi contains all the basis
  functions and C has c_j's as columns.

  The second optimization problem that can be solved is the following

  min_i min t_i

  subject to

            t0 <= t_i <= t0 + eta
            A_i^T p(t_i) = b_i
            A_j^T p(t_i) <= b_i         j neq i



  We use interior-point methods to solve the non-linear optimization program
  using COIN-OR ipopt and the ifopt interface for Eigen + ipopt.

*/

#ifndef NLP_HPOLYORACLES_HPP
#define NLP_HPOLYORACLES_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>

#define POLYTOL 1e-7

using namespace ifopt;

// Define the variable t we use in the optimization
template <typename VT, typename NT>
class HPolyOracleVariables : public VariableSet {
public:
  NT t, tb, eta;

  HPolyOracleVariables(NT t_prev, NT tb_, NT eta_=-1): VariableSet(1, "t"), t(t_prev), tb(tb_), eta(eta_) {};

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
    NT tu = eta > 0 ? NT(tb + eta) : NT(inf);

    bounds.at(0) = Bounds(tb, tu);
    return bounds;
  };

};

// Define the cost function f(t) = t (ipopt takes minimization so it is -t)
template <typename VT, typename NT>
class HPolyOracleCost : public CostTerm {
public:
  std::string method;

  HPolyOracleCost(std::string method_) : CostTerm("h_poly_cost"), method(method_) {};

  NT GetCost() const override {
    VectorXd T = GetVariables()->GetComponent("t")->GetValues();
    if (method == "max_pos") return - T(0);
    else if (method == "min_pos") return T(0);
  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac) const override {
    if (var_set == "t") {
      if (method == "max_pos") jac.coeffRef(0, 0) = (NT) (-1.0);
      else if (method == "min_pos")  jac.coeffRef(0, 0) = (NT) (1.0);
    }
  }

};

// Define the feasibility constraint A p(t) <= b which we translate
// to (A * C) * Phi <= b
template <typename MT, typename VT, typename NT, class bfunc>
class HPolyOracleFeasibility : public ConstraintSet {
public:
  MT &C;
  VT &b;
  bfunc phi, grad_phi;
  NT t0;
  int m, M;
  std::string method;
  int index;

  HPolyOracleFeasibility(MT &C_, VT &b_, NT t0_, bfunc basis, bfunc basis_grad, std::string method_, int i) :
    C(C_), b(b_), t0(t0_), phi(basis), grad_phi(basis_grad), ConstraintSet(C_.rows(), "h_poly_feasibility"),
    method(method_), index(i) {
      m = C_.rows();
      M = C_.cols();
    };

  // Define bounds for feasibility
  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    for (int i = 0; i < m; i++) {
      if (method == "min_pos" && i == index) {
        bounds.at(i) = Bounds(b(i), b(i));
      }
      else {
        bounds.at(i) = Bounds((NT) (-inf), b(i));
      }
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

  // Calculate jacobian matrix of constraints
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

// Helper function that calls the optimization problem (called from hpolytope.h)
template <typename MT, typename VT, typename Point, typename NT, class bfunc>
std::tuple<NT, Point, int> curve_intersect_hpoly_ipopt_helper(NT t_prev, NT t0, NT eta, MT &A, VT &b, std::vector<Point> &coeffs, bfunc phi, bfunc grad_phi, std::string solution="max_pos")
{


  MT C, C_tmp;
  C_tmp.resize(coeffs[0].dimension(), coeffs.size());


  for (int i = 0; i < coeffs.size(); i++) {
    C_tmp.col(i) = coeffs[i].getCoefficients();
  }

  // C_tmp: dimension x num_coeffs
  // A: constraints x dimension
  // C: constraints x num_coeffs
  C = A * C_tmp;

  // Initialize COIN-OR ipopt solver
  IpoptSolver ipopt;
  ipopt.SetOption("linear_solver", "mumps");
  ipopt.SetOption("jacobian_approximation", "exact");
  ipopt.SetOption("tol", POLYTOL);
  ipopt.SetOption("acceptable_tol", 100 * POLYTOL);
  ipopt.SetOption("max_iter", 1000000);

  ipopt.SetOption("print_level", 0);
  ipopt.SetOption("sb", "yes");

  NT t, t_tmp;

  if (solution == "max_pos") {

      Problem nlp;
      std::shared_ptr<HPolyOracleVariables<VT, NT>> hpolyoraclevariables (new HPolyOracleVariables<VT, NT>(t_prev, t0, eta));
      std::shared_ptr<HPolyOracleCost<VT, NT>> hpolyoraclecost (new HPolyOracleCost<VT, NT>(solution));
      std::shared_ptr<HPolyOracleFeasibility<MT, VT, NT, bfunc>> hpolyoraclefeasibility (new HPolyOracleFeasibility<MT, VT, NT, bfunc>(C, b, t0, phi, grad_phi, "max_pos", 0));

      nlp.AddVariableSet  (hpolyoraclevariables);
      nlp.AddCostSet      (hpolyoraclecost);
      nlp.AddConstraintSet(hpolyoraclefeasibility);

      ipopt.Solve(nlp);

      t = nlp.GetOptVariables()->GetValues()(0);

  }

  if (solution == "min_pos") {
    int m = A.rows();

    t = eta > 0 ? t0 + eta : std::numeric_limits<NT>::max();

    for (int i = 0; i < m; i++) {

      Problem nlp;
      std::shared_ptr<HPolyOracleVariables<VT, NT>> hpolyoraclevariables (new HPolyOracleVariables<VT, NT>(t_prev, t0, eta));
      std::shared_ptr<HPolyOracleCost<VT, NT>> hpolyoraclecost (new HPolyOracleCost<VT, NT>(solution));
      std::shared_ptr<HPolyOracleFeasibility<MT, VT, NT, bfunc>> hpolyoraclefeasibility (new HPolyOracleFeasibility<MT, VT, NT, bfunc>(C, b, t0, phi, grad_phi, "min_pos", i));

      nlp.AddVariableSet  (hpolyoraclevariables);
      nlp.AddCostSet      (hpolyoraclecost);
      nlp.AddConstraintSet(hpolyoraclefeasibility);

      t_tmp = nlp.GetOptVariables()->GetValues()(0);

      std::cout << "t is " << t_tmp << std::endl;

      if (t_tmp < t && t_tmp > t0) t = t_tmp;


    }
  }

  Point p(coeffs[0].dimension());

  for (unsigned int i = 0; i < coeffs.size(); i++) {
    p += phi(t, t0, i, coeffs.size()) * coeffs[i];
  }

  const NT* b_data = b.data();

  int f_min = -1;
  NT dist_min = std::numeric_limits<NT>::max();

  for (int i = 0; i < A.rows(); i++) {
    if (*b_data == 0 && std::abs(*b_data - A.row(i) * p.getCoefficients()) < dist_min) {
      f_min = i;
      dist_min = std::abs(*b_data - A.row(i) * p.getCoefficients());
    }
    else if (*b_data != 0 && std::abs(*b_data - A.row(i) * p.getCoefficients()) / *b_data < dist_min) {
      f_min = i;
      dist_min = std::abs(*b_data - A.row(i) * p.getCoefficients()) / *b_data;
    }
    b_data++;
  }

  return std::make_tuple(t, p, f_min);
}

#endif
