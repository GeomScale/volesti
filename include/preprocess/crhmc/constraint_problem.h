// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONSTRAINT_PROBLEM_H
#define CONSTRAINT_PROBLEM_H
#include "Eigen/Eigen"
/*Input structure: With this the user can define a polytope sampling problem*/
template <typename MatrixType, typename Point>
class constraint_problem {
public:
  using Type = typename Point::FT;
  using PointType = Point;
  using VT = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
  using MT = MatrixType;
private:
  unsigned int num_vars; // num_vars of the original problem
  MatrixType Aeq;         // Matrix of coefficients for the equality constraints
  VT beq;                 // Right hand side of the equality constraints
  MatrixType Aineq; // Matrix of coefficients for the inequality constraints
  VT bineq;         // Right hand side of the inequality constraints
  VT lb;            // lb on the output coordinates preset to -1e9
  VT ub;            // ub on the output coordinates preset to +1e9
  Type inf = 1e9;
public:
  /*Constructors for different input instances*/
  constraint_problem(const int dim, MT const &Aeq_, VT const &beq_, MT const &Aineq_,
              VT const &bineq_, VT const &lb_, VT const &ub_)
      : num_vars(dim), Aeq(Aeq_), beq(beq_), Aineq(Aineq_), bineq(bineq_),
        lb(lb_), ub(ub_) {
        }

  constraint_problem(const int dim) {
    num_vars = dim;
    init(num_vars);
  }


  void init(int num_vars) {
    Aineq.resize(0, num_vars);
    Aeq.resize(0, num_vars);
    bineq.resize(0, 1);
    beq.resize(0, 1);
    lb = -VT::Ones(num_vars) * inf;
    ub = VT::Ones(num_vars) * inf;
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
  std::pair<MT,VT> get_equations(){
    return std::make_pair(Aeq,beq);
  }
  std::pair<MT,VT> get_inequalities(){
    return std::make_pair(Aineq,bineq);
  }
  std::pair<VT,VT> get_bounds(){
    return std::make_pair(lb,ub);
  }
  unsigned int dimension(){
    return num_vars;
  }

};

#endif
