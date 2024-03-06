// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_ORACLE_FUNCTORS_RCPP_HPP
#define ODE_SOLVERS_ORACLE_FUNCTORS_RCPP_HPP

enum ode_solvers {
  no_solver,
  leapfrog,
  euler,
  runge_kutta,
  richardson,
  collocation,
  implicit_midpoint
};

// Holds Oracle Functor that wraps an R function via RCpp
// The R function is provided as an Rcpp::Function object
// The functor uses Rcpp::as and Rcpp::wrap to do the conversion,
// call the oracle, and convert the results back to C++
struct RcppFunctor {

  template <
      typename NT
  >
  struct parameters {
    NT L; // Lipschitz constant of gradient
    NT m; // Strong-convexity parameter
    NT eta; // Step-size (if defined by user)
    NT kappa; // Condition number
    unsigned int order; // Order of ODE functor

    parameters(
      NT L_,
      NT m_,
      NT eta_,
      unsigned int order_=2
    ) :
      L(L_),
      m(m_),
      eta(eta_),
      kappa(L_ / m_),
      order(order_)
    {}
  };

  // Log-probability gradient functor
  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef typename Point::Coeff VT;
    typedef std::vector<Point> pts;

    parameters<NT> params;
    Rcpp::Function neg_grad_f; // Negative gradient as an Rcpp::Function
    bool negate;

    GradientFunctor(
      parameters<NT> params_,
      Rcpp::Function neg_grad_f_,
      bool negate_=true):
      params(params_),
      neg_grad_f(neg_grad_f_),
      negate(negate_)
    {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        // Convert point to Rcpp::NumericMatrix

        VT y = Rcpp::as<VT>(neg_grad_f(Rcpp::wrap(xs[0].getCoefficients())));

        Point z(y);

        if (negate) z = (-1.0) * z;

        // Return result as Point
        return z;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

    Point operator() (Point const& x) const {
      VT y = Rcpp::as<VT>(neg_grad_f(Rcpp::wrap(x.getCoefficients())));

      Point z(y);

      if (negate) z = (-1.0) * z;

      // Return result as Point
      return z;
    }

  };

  // Negative log-probability functor
  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;
    typedef typename Point::Coeff VT;

    parameters<NT> params;
    Rcpp::Function negative_logprob;

    FunctionFunctor(
      parameters<NT> params_,
      Rcpp::Function negative_logprob_) :
      params(params_),
      negative_logprob(negative_logprob_)
    {};

    NT operator() (Point const& x) const {
      return Rcpp::as<NT>(negative_logprob(Rcpp::wrap(x.getCoefficients())));
    }

  };

  // Log-probability hessian functor
  template
  <
    typename Point
  >
  struct HessianFunctor {
    typedef typename Point::FT NT;
    typedef typename Point::Coeff VT;

    parameters<NT> params;
    Rcpp::Function hessian; // Negative hessian as an Rcpp::Function

    HessianFunctor(
      parameters<NT> params_,
      Rcpp::Function hessian_) :
      params(params_),
      hessian(hessian_)
    {};

    Point operator() (Point const& x) const {
      VT y= Rcpp::as<VT>(hessian(Rcpp::wrap(x.getCoefficients())));
      return Point(y);
    }

  };
};

#endif
