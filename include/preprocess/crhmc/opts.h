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

#ifndef OPTS_H
#define OPTS_H

/// @brief Crhmc options
/// @tparam Type Numer type
template <typename Type> class opts {
public:
  /*Preprocess options*/
  const int ipmMaxIter = 200; //Maximum number of iterations for finding the analytic and lewis center
  const Type ipmDistanceTol = 1e-8;
  const Type ipmDualTol = 1e-12;
  int maxNZ = 30;
  Type max_coord = 1e9;
  bool EnableReordering = true;
  const int min_convergence_steps=8;

  /*ODE options*/
  const Type implicitTol = 1e-5;
  const int maxODEStep = 30;
  Type initialStep = 0.2;
  Type convergence_bound = 1e16;

  /*PackedCS Solver Options*/
  Type solver_accuracy_threshold=1e-2;
  int simdLen=1;

  /*Sampler options*/
  bool DynamicWeight = true; //Enable the use of dynamic weights for each variable when sampling
  bool DynamicStepSize = true;  // Enable adaptive step size that avoids low acceptance probability
  bool DynamicRegularizer = true; //Enable the addition of a regularization term
  Type regularization_factor=1e-20;
  /*Dynamic step choices*/
  Type warmUpStep = 10;
  int maxConsecutiveBadStep = 10;
  Type targetODEStep = 10;
  Type shrinkFactor = 1.1;
  Type minStepSize = 1e-5;
  Type effectiveStepSize = 1;

  opts() {}
  void operator=(const opts &rhs) {
    EnableReordering = rhs.EnableReordering;
    maxNZ = rhs.maxNZ;
  }
};
#endif
