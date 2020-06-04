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

using namespace ifopt;

template <typename VT, typename Point, typename NT>
class HPolyOracleVariables : VariableSet {
public:
  NT t, tb;

  HPolyOracleVariables(NT t0, NT tb_=NT(0)): VariableSet(1, "t"), t(t0), tb(tb_) {};

  void SetVariables(const VT& T) override {
    t = T(0);
  }

  VT GetValues() const override {
    VT T;
    VT.resize(1);
    T(0) = t;
    return T;
  }

  NT GetValue() const override {
    return t;
  }

  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    bounds.at(0) = Bounds(tb, (NT) inf);
    return bounds;
  };

};

template <typename VT, typename Point, typename NT>
class HPolyOracleCost : public CostTerm {
  HPolyOracleCost() : CostTerm("h_poly_cost") {};

  NT GetCost() const override {
    GetVariables()->GetComponent("t")->GetValue();
  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac) const override {
    if (var_set == "t") jac.coeffRef(0, 0) = (NT) (1.0);
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
      bounds.at(i) = b(i);
    }
    return bounds;
  }

  VT GetValues() const override {
    NT t = GetVariables()->GetComponent("t")->GetValue();
    VT phis;
    phis.resize(M);

    for (int i = 0; i < M; i++) {
      phis(i) = (NT) (phi(t, t0, i, M));
    }
    return C * phis;
  }

  void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override {

    if (var_set == "t") {
      NT t = GetVariables()->GetComponent("t")->GetValue();
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
std::tuple<NT, Point, int> curve_intersect_ipopt(NT t0, const MT &A, const VT &b, std::vector<Point> &coeffs, bfunc phi, bfunc grad_phi)
{

  Problem nlp;

  nlp.AddVariableSet  (std::make_shared<HPolyOracleVariables>());
  nlp.AddCostSet      (std::make_shared<HPolyOracleCost>());

  IpoptSolver solver;
  ipopt.SetOption("linear_solver", "mumps");
  ipopt.SetOption("jacobian_approximation", "exact");

  ipopt.Solve(nlp);

  NT t = nlp.GetOptVariables()->GetValue();

  Point p(coeffs[0].dimension());

  for (unsigned int i = 0; i < coeffs.size(); i++) {
    p += phi(t, t0, i, coeffs.size()) * coeffs[i];
  }

  return std::make_tuple(t, p, -1);

}


#endif
