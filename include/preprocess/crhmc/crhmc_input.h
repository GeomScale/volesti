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

template <typename MatrixType, typename Type> class crhmc_input {
  using VT = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
  const Type inf = 1e9;

public:
  MatrixType Aineq;
  VT bineq;
  MatrixType Aeq;
  VT beq;
  opts<Type> options;
  VT lb;
  VT ub;
  crhmc_input(int dimension) {
    Aineq.resize(0, dimension);
    Aeq.resize(0, dimension);
    bineq.resize(0, 1);
    beq.resize(0, 1);
    lb = -VT::Ones(dimension) * inf;
    ub = VT::Ones(dimension) * inf;
  }
};
#endif
