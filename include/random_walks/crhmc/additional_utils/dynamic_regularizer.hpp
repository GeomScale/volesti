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
#ifndef DYNAMIC_REGULARIZER_HPP
#define DYNAMIC_REGULARIZER_HPP
template <typename NT> class dynamic_regularizer {
  NT bound = 1;
  dynamic_regularizer(Hamiltonian &ham) : { ham.barrier.extraHessian = 1; }
  void step(Hamiltonian &ham) {
    bound = (ham.x.getCoefficients().maxCoeff()), bound);
    if (!ham.freezed) {
    }
  }
};
#endif
