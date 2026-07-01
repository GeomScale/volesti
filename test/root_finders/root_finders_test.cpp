// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DISABLE_NLP_ORACLES

#include <cstdlib>
#include <unistd.h>
#include <pthread.h>
#include <mps/mps.h>
#include <gmp.h>
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include <string>
#include <typeinfo>

#include "Eigen/Eigen"
#include "doctest.h"

#include "root_finders.hpp"

template<typename NT>
void test_newton_raphson() {
  typedef std::function<NT(NT)> func;

  func f = [](NT t) {
    return t * t - 2 * t + 1;
  };

  func grad_f = [](NT t) {
    return 2 * t - 2;
  };

  NT t0 = 5.0;

  NT t = newton_raphson<NT, func>(t0, f, grad_f, 0.001).first;

  std::cout << t << std::endl;

  CHECK(std::abs(t - 1) < 0.001);

}

template<typename NT>
void test_mpsolve() {
  NT S = NT(5);
  NT P = NT(6);

  CHECK(S >= NT(0));
  CHECK(P > NT(0));

  std::vector<NT> coeffs{P, -S, NT(1)};

  std::vector<std::pair<NT, NT>> results = mpsolve<NT>(coeffs, true);

  NT x1 = results[0].first;
  NT x2 = results[1].first;

  NT S_ = x1 + x2;
  NT P_ = x1 * x2;

  CHECK(std::abs(S_ - S) / S < 0.001);

  CHECK(std::abs(P_ - P) / P < 0.001);

}

template<typename NT>
void call_test_root_finders() {
  std::cout << "--- Testing Newton-Raphson" << std::endl;
  test_newton_raphson<NT>();

  std::cout << "--- Test mpsolve" << std::endl;
  test_mpsolve<NT>();

}

TEST_CASE("root_finders") {
  call_test_root_finders<double>();
}

#endif
