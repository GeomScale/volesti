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
template <typename Type> class opts {
public:
  /*Preprocess options*/
  const int ipmMaxIter = 200;
  const Type ipmDistanceTol = 1e-8;
  const Type ipmDualTol = 1e-12;
  int maxNZ = 30;
  Type max_coord = 1e7;
  bool EnableReordering = true;

  /*ODE options*/
  const Type implicitTol = 1e-5;
  const int maxODEStep = 30;
  Type initialStep = 0.2;

  /*Sampler options*/
  bool DynamicWeight = false;
  bool DynamicStepSize = false;
  bool DynamicRegularizer = false;

  /*Dynamic step choices*/
  Type warmUpStep = 10;
  Type maxConsecutiveBadStep = 10;
  Type targetODEStep = 10;
  Type shrinkFactor = 1.1;
  Type minStepSize = 0.001;
  Type effectiveStepSize = 1;

  opts() {}
  void operator=(const opts &rhs) {
    EnableReordering = rhs.EnableReordering;
    maxNZ = rhs.maxNZ;
  }
};
#endif
