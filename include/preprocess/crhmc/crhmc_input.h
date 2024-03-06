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
#ifndef CRHMC_INPUT_H
#define CRHMC_INPUT_H
#include "Eigen/Eigen"
#include "opts.h"
#include "convex_bodies/hpolytope.h"
#include "preprocess/crhmc/constraint_problem.h"
/*0 funciton handles are given as a reference in case the user gives no
function. Then the uniform function is implied*/
template <typename Point>
struct ZeroFunctor
{
  using Type = typename Point::FT;
  Point operator()(Point const &x) const { return Point(x.dimension()); }
  struct parameters {
    Type L=1;
    Type eta=1;
  };
  struct parameters params;
};
template <typename Point>
struct ZeroScalarFunctor
{
  using Type = typename Point::FT;
  Type operator()(Point const &x) const { return 0; }
};

///
/// Input structure: With this the user can define the input for a crhmc polytope sampling problem
template <typename MatrixType, typename Point,
          typename func = ZeroScalarFunctor<Point>,
          typename grad = ZeroFunctor<Point>,
          typename hess = ZeroFunctor<Point>>
class crhmc_input
{
  using Type = typename Point::FT;
  using VT = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
  ZeroFunctor<Point> zerof;
  ZeroScalarFunctor<Point> zerosf;

public:
  using MT = MatrixType;
  using Func = func;
  using Grad = grad;
  using Hess = hess;
  using point= Point;
  MatrixType Aineq;                       // Matrix of coefficients for the inequality constraints
  VT bineq;                               // Right hand side of the inequality constraints
  MatrixType Aeq;                         // Matrix of coefficients for the equality constraints
  VT beq;                                 // Right hand side of the equality constraints
  opts<Type> options;                     // structure of the parameters of the problem
  VT lb;                                  // lb on the output coordinates preset to -1e7
  VT ub;                                  // ub on the output coordinates preset to +1e7
  func &f;                                // Negative log density function handle
  grad &df;                               // Negative log density gradient function handle
  hess &ddf;                              // Negative log density hessian function handle
  bool fZero;                             // whether f is completely zero
  bool fHandle;                           // whether f is handle or not
  bool dfHandle;                          // whether df is handle or not
  bool ddfHandle;                         // whether ddf is handle or not
  unsigned int dimension;                 // dimension of the original problem
  const Type inf = options.max_coord + 1; // helper for barrier handling
  /*Constructors for different input instances*/
  crhmc_input(int dim, func &function, grad &g, hess &h)
      : f(function), df(g), ddf(h)
  { dimension=dim;
    fZero = false;
    fHandle = true;
    dfHandle = true;
    ddfHandle = true;
    init(dimension);
  }
  crhmc_input(int dim, func &function)
      : f(function), df(zerof), ddf(zerof)
  { dimension=dim;
    fZero = false;
    fHandle = true;
    dfHandle = false;
    ddfHandle = false;
    init(dimension);
  }
  crhmc_input(int dim, func &function, grad &g)
      : f(function), df(g), ddf(zerof)
  { dimension=dim;
    fZero = false;
    fHandle = true;
    dfHandle = true;
    ddfHandle = false;
    init(dimension);
  }
  crhmc_input(int dim) : f(zerosf), df(zerof), ddf(zerof)
  { dimension=dim;
    fZero = true;
    fHandle = false;
    dfHandle = false;
    ddfHandle = false;
    init(dimension);
  }

  void init(int dimension)
  {
    Aineq.resize(0, dimension);
    Aeq.resize(0, dimension);
    bineq.resize(0, 1);
    beq.resize(0, 1);
    lb = -VT::Ones(dimension) * inf;
    ub = VT::Ones(dimension) * inf;
  }
};
#include <type_traits>

template <
    typename Input, typename Polytope, typename Func, typename Grad, typename Hess,
    typename std::enable_if<std::is_same<
        Polytope, HPolytope<typename Input::point>>::value>::type * = nullptr>
inline Input convert2crhmc_input(Polytope &P, Func &f, Grad &g, Hess &h) {
  int dimension = P.dimension();
  Input input = Input(dimension, f, g, h);
  if (std::is_same<
      Hess, ZeroFunctor<typename Input::point>>::value){
        input.ddfHandle=false;
  }
  input.Aineq = P.get_mat();
  input.bineq = P.get_vec();
  return input;
}

template <typename Input, typename Polytope, typename Func, typename Grad, typename Hess,
          typename std::enable_if<std::is_same<
              Polytope, constraint_problem<typename Input::MT,
                                           typename Input::point>>::value>::type
              * = nullptr>
inline Input convert2crhmc_input(Polytope &P, Func &f, Grad &g, Hess &h) {
  int dimension = P.dimension();
  Input input = Input(dimension, f, g, h);
  if (std::is_same<
      Hess, ZeroFunctor<typename Input::point>>::value){
        input.ddfHandle=false;
  }
  std::tie(input.Aineq, input.bineq) = P.get_inequalities();
  std::tie(input.Aeq, input.beq) = P.get_equations();
  std::tie(input.lb, input.ub) = P.get_bounds();
  return input;
}
template <typename Input, typename Polytope, typename Func, typename Grad, typename Hess,
          typename std::enable_if<
              !std::is_same<Polytope,
                            constraint_problem<typename Input::MT,
                                               typename Input::point>>::value &&
              !std::is_same<Polytope, HPolytope<typename Input::point>>::value>::
              type * = nullptr>
inline Input convert2crhmc_input(Polytope &P, Func &f, Grad &g, Hess &h) {
  /*CRHMC works only for H-polytopes and constraint_problems for now*/
  int dimension = 0;
  Input input = Input(dimension, f, g, h);
  return input;
}
#endif
