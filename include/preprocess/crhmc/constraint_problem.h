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
#ifndef CONSTRAINT_PROBLEM_H
#define CONSTRAINT_PROBLEM_H
#include "Eigen/Eigen"
/*Input structure: With this the user can define a polytope sampling problem*/
template <typename MatrixType, typename Point>
class constraint_problem {
public:
  using Type = typename Point::FT;
  using VT = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
  using MT = MatrixType;

  using point = Point;
  unsigned int dimension; // dimension of the original problem
  MatrixType Aeq;         // Matrix of coefficients for the equality constraints
  VT beq;                 // Right hand side of the equality constraints
  MatrixType Aineq; // Matrix of coefficients for the inequality constraints
  VT bineq;         // Right hand side of the inequality constraints
  VT lb;            // lb on the output coordinates preset to -1e7
  VT ub;            // ub on the output coordinates preset to +1e7
  Type inf = 1e7 + 1;
  /*Constructors for different input instances*/
  constraint_problem(const int dim, MT const &Aeq_, VT const &beq_, MT const &Aineq_,
              MT const &bineq_, VT const &lb_, VT const &ub_)
      : dimension(dim), Aeq(Aeq_), beq(beq_), Aineq(Aineq_), bineq(bineq_),
        lb(lb_), ub(ub_) {
        }

  constraint_problem(const int dim) {
    dimension = dim;
    init(dimension);
  }


  void init(int dimension) {
    Aineq.resize(0, dimension);
    Aeq.resize(0, dimension);
    bineq.resize(0, 1);
    beq.resize(0, 1);
    lb = -VT::Ones(dimension) * inf;
    ub = VT::Ones(dimension) * inf;
  }
  void set_equality_constraints(MT const &Aeq_, VT const &beq_){
    Aeq = Aeq_;
    beq = beq_;
  }
  void set_inequality_constraints(MT const &Aineq_, VT const &bineq_){
    Aineq = Aineq_;
    bineq = bineq_;
  }
  void set_bounds(VT const &lb_, VT const &ub_){
    lb = lb_;
    ub = ub_;
  }
};

#endif
