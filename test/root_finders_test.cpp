/*
GeomScale Project

Copyright (c) 2020
  Vissarion Fisikopoulos
  Apostolos Chalkis
  Elias Tsigaridas
  Marios Papachristou

Contributed and/or modified by Marios Papachristou,
as part of Google Summer of Code 2020 program.

VolEsti is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

VolEsti is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

See the file COPYING.LESSER for the text of the GNU Lesser General
Public License.  If you did not receive this file along with HeaDDaCHe,
see <http://www.gnu.org/licenses/>.
*/

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
