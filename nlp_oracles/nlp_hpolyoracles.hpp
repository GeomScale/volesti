// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

/*  This file contains the implementation of boundary oracles between an
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
#include <limits>
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>

#include "root_finders/root_finders.hpp"
#define MAX_NR_TRIES 10000

#define POLYTOL 1e-7

#ifndef isnan
  using std::isnan;
#endif

#ifndef isinf
  using std::isinf;
#endif


using namespace ifopt;

/// Define the variable t we use in the optimization
/// \tparam VT Vector Type
/// \tparam NT Numeric Type
template <typename VT, typename NT>
class HPolyOracleVariables : public VariableSet {
public:
  NT t, tb, eta;

  HPolyOracleVariables(NT t_prev, NT tb_, NT eta_=-1):
    VariableSet(1, "t"), t(t_prev), tb(tb_), eta(eta_) {};

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

/// Define the cost function f(t) = t (ipopt takes minimization so it is -t)
/// \tparam VT Vector Type
/// \tparam NT Numeric Type
template <typename VT, typename NT>
class HPolyOracleCost : public CostTerm {
public:
  std::string method;

  HPolyOracleCost(std::string method_) :
    CostTerm("h_poly_cost"), method(method_) {};

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

/// Define the feasibility constraint A p(t) <= b which we translate
/// to (A * C) * Phi <= b
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
/// \tparam NT Numeric Type
/// \tparam bfunc feasibility constraint type
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
  int ignore_facet;

  HPolyOracleFeasibility(
    MT &C_,
    VT &b_,
    NT t0_,
    bfunc basis,
    bfunc basis_grad,
    std::string method_,
    int i,
    int ignore_facet_=-1) :
    C(C_), b(b_), t0(t0_), phi(basis), grad_phi(basis_grad),
    ConstraintSet(C_.rows(), "h_poly_feasibility"),
    method(method_), index(i), ignore_facet(ignore_facet_) {
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

      if (i == ignore_facet) {
        bounds.at(i) = Bounds((NT) (-inf), (NT) (inf));
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
          if (isinf(temp)) continue;
          else jac_block.coeffRef(i, 0) += temp;
        }
      }
    }
  }


};

/// Oracle that uses the COIN-OR ipopt solver
/// \tparam Polytope Polytope Type
/// \tparam bfunc feasibility constraint type
template <typename Polytope, class bfunc>
struct IpoptHPolyoracle {
  typedef typename Polytope::MT MT;
  typedef typename Polytope::VT VT;
  typedef typename Polytope::NT NT;
  typedef typename Polytope::PointType Point;


  std::tuple<NT, Point, int> apply(
    NT t_prev,
    NT t0,
    NT eta,
    MT &A,
    VT &b,
    Polytope &P,
    std::vector<Point> &coeffs,
    bfunc phi,
    bfunc grad_phi,
    int ignore_facet=-1,
    std::string solution="max_pos")
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
        std::shared_ptr<HPolyOracleVariables<VT, NT>>
          hpolyoraclevariables (new HPolyOracleVariables<VT, NT>(t_prev, t0, eta));
        std::shared_ptr<HPolyOracleCost<VT, NT>>
          hpolyoraclecost (new HPolyOracleCost<VT, NT>(solution));
        std::shared_ptr<HPolyOracleFeasibility<MT, VT, NT, bfunc>>
          hpolyoraclefeasibility (new HPolyOracleFeasibility<MT, VT, NT, bfunc>
            (C, b, t0, phi, grad_phi, "max_pos", 0, ignore_facet));

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

        if (i == ignore_facet) continue;

        Problem nlp;
        std::shared_ptr<HPolyOracleVariables<VT, NT>>
          hpolyoraclevariables (new HPolyOracleVariables<VT, NT>(t_prev, t0, eta));
        std::shared_ptr<HPolyOracleCost<VT, NT>>
          hpolyoraclecost (new HPolyOracleCost<VT, NT>(solution));
        std::shared_ptr<HPolyOracleFeasibility<MT, VT, NT, bfunc>>
          hpolyoraclefeasibility (new HPolyOracleFeasibility<MT, VT, NT, bfunc>
            (C, b, t0, phi, grad_phi, "min_pos", i));

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
  };

};

/// Compute intersection of H-polytope P := Ax <= b
/// with polynomial curve p(t) = sum a_j (t - t0)^j
/// Uses the MPsolve library
/// \tparam Polytope Polytope Type
/// \tparam bfunc feasibility constraint type
template <typename Polytope, class bfunc>
struct MPSolveHPolyoracle {

  typedef typename Polytope::MT MT;
  typedef typename Polytope::VT VT;
  typedef typename Polytope::NT NT;
  typedef typename Polytope::PointType Point;

  std::tuple<NT, Point, int> apply(
    NT t_prev,
    NT t0,
    NT eta,
    MT &A,
    VT &b,
    Polytope &P,
    std::vector<Point> &coeffs,
    bfunc phi,
    bfunc grad_phi,
    int ignore_facet=-1,
    bool positive_real=true)
  {
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();

    NT tu = eta > 0 ? t0 + eta : NT(maxNT);
    NT t = tu;
    Point dummy(coeffs[0].dimension());

    for (unsigned int j = 0; j < coeffs.size(); j++) {
      dummy = dummy + pow(tu - t0, NT(j)) * coeffs[j];
    }

    std::tuple<NT, Point, int> result = std::make_tuple(tu, dummy, -1);

    int m = A.rows();

    // Keeps constants A_i^T C_j
    std::vector<NT> Z(coeffs.size(), NT(0));

    // std::vector<std::pair<NT, NT>> solutions;

    // Iterate over all hyperplanes
    for (int i = 0; i < m; i++) {

      if (i == ignore_facet) continue;

      for (unsigned int j = 0; j < coeffs.size(); j++) {
        Z[j] = A.row(i) * coeffs[j].getCoefficients();

        #ifdef VOLESTI_DEBUG
        std::cout << "Z [ " << j << " ] = " << Z[j] << std::endl;
        #endif
      }

      // Find point projection on m-th hyperplane
      Z[0] -= b(i);

      std::vector<std::pair<NT, NT>> solutions = mpsolve<NT>(Z, positive_real);

      for(std::pair<NT, NT> sol: solutions) {

        #ifdef VOLESTI_DEBUG
        std::cout << "Facet: " << i << " Candidate root is " << sol.first + t0 << std::endl;
        #endif

        // Check if solution is in the desired range [t0, t0 + eta] and if it is the current minimum
        if (t0 + sol.first <= tu && t0 + sol.first < std::get<0>(result)) {
          t = t0 + sol.first;

          // Calculate point from this root
          Point p = Point(coeffs[0].dimension());

          for (unsigned int j = 0; j < coeffs.size(); j++) {
            p += pow(t - t0, NT(j)) * coeffs[j] ;
          }

          #ifdef VOLESTI_DEBUG
          std::cout << "Calculcated point is " << std::endl << p.getCoefficients() << std::endl;
          #endif

          // Check if point satisfies Ax <= b up to some tolerance and change current solution
          if (P.is_in(p, 1e-6)) {
            result =  std::make_tuple(t, p, i);
          }
        }
      }

    }

    return result;
  };

};

/// Compute intersection of H-polytope P := Ax <= b
/// with curve p(t) = sum a_j phi_j(t) where phi_j are basis
/// functions (e.g. polynomials)
/// Uses Newton-Raphson to solve the transcendental equation
/// \tparam Polytope Polytope Type
/// \tparam bfunc feasibility constraint type
template <typename Polytope, class bfunc>
struct NewtonRaphsonHPolyoracle {
  typedef typename Polytope::MT MT;
  typedef typename Polytope::VT VT;
  typedef typename Polytope::NT NT;
  typedef typename Polytope::PointType Point;

  std::tuple<NT, Point, int> apply(
    NT t_prev,
    NT t0,
    NT eta,
    MT &A,
    VT &b,
    Polytope &P,
    std::vector<Point> &coeffs,
    bfunc phi,
    bfunc grad_phi,
    int ignore_facet=-1)
  {

    // Keep results in a vector (in case of multiple roots)
    // The problem has O(m * len(coeffs)) solutions if phi's are polys
    // due to the Fundamental Theorem of Algebra
    // Some roots may be common for more than one hyperplanes
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();

    // Root
    NT t = t_prev;
    NT tu = eta > 0 ? t0 + eta : NT(maxNT);

    Point dummy(coeffs[0].dimension());

    for (unsigned int j = 0; j < coeffs.size(); j++) {
      dummy += coeffs[j] * phi(tu, t0, j, coeffs.size());
    }


    std::tuple<NT, Point, int> result = std::make_tuple(tu, dummy, -1);

    // Helper variables for Newton-Raphson
    NT dot_u, num, den, den_tmp;

    // Regularization for NR (e.g. at critical points where grad = 0)
    NT reg = (NT) 1e-7;
    VT Z;
    int m = A.rows();

    // Keeps constants A_i^T C_j
    Z.resize(coeffs.size());

    // Iterate over all hyperplanes
    for (int i = 0; i < m; i++) {

      if (i == ignore_facet) continue;

      // Calculate constants
      start_iter: t_prev = t0 + reg;


      for (unsigned int j = 0; j < coeffs.size(); j++) {
        Z(j) = A.row(i) * coeffs[j].getCoefficients();
      }


      for (int tries = 0; tries < MAX_NR_TRIES; tries++) {

        num = - b(i);
        den = (NT) 0;

        // Calculate numerator f(t) and denominator f'(t)
        for (int j = 0; j < coeffs.size(); j++) {
          num += Z(j) * phi(t_prev, t0, j, coeffs.size());

          // Avoid ill-posed derivative (e.g. 0^{-1})
          if (j > 0) den += Z(j) * grad_phi(t_prev, t0, j, coeffs.size());
        }

        // Regularize denominator if near 0
        if (std::abs(den) < 10 * reg) den += reg;

        // Newton-Raphson Iteration t = t_old - f(t) / f'(t)
        t = t_prev - num / den;

        if (t < 0 && t_prev < 0) continue;

        if (std::abs(t - t_prev) < 1e-6 && t > t0) {
          // Add root (as t) and facet

          Point p = Point(coeffs[0].dimension());

          for (unsigned int j = 0; j < coeffs.size(); j++) {
            p += coeffs[j] * phi(t, t0, j, coeffs.size());
          }

          if (P.is_in(p) && t < std::get<0>(result))
            result =  std::make_tuple(t, p, i);

        }

        t_prev = t;

      }

    }

    return result;

  };

};


#endif
