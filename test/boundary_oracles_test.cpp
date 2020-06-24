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

#include "Eigen/Eigen"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"

#include "known_polytope_generators.h"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

#include "exact_vols.h"
#include "generators/known_polytope_generators.h"

#include <string>
#include <typeinfo>
#include <chrono>
#include "doctest.h"

template <typename NT, class Point, class bfunc>
void test_h_poly_oracles(std::vector<Point> coeffs, bfunc phi, bfunc grad_phi, NT t_des, int facet_des) {
  typedef boost::mt19937    RNGType;
  typedef HPolytope<Point> Hpolytope;
  typedef std::tuple<NT, Point, int> result;
  Hpolytope P;
  NT tol = 1e-4;

  P = gen_cube<Hpolytope>(2, false);

  result res = P.curve_intersect_newton_raphson(0.01, 0, -1, coeffs, phi, grad_phi);
  NT t = std::get<0>(res);
  int facet = std::get<2>(res);

  // CHECK(facet == facet_des);
  CHECK(std::abs(t - t_des) / t_des < tol);

  result res2 = P.curve_intersect_ipopt(0.0, 0, -1, coeffs, phi, grad_phi);

  t = std::get<0>(res2);
  facet = std::get<2>(res2);
  CHECK(std::abs(t - t_des) / t_des < tol);

  res2 = P.curve_intersect_mpsolve(0, 0, -1, coeffs);

  t = std::get<0>(res2);
  // std::cout << "t is " << t << std::endl;
  facet = std::get<2>(res2);
  CHECK(std::abs(t - t_des) / t_des < tol);


}

template <typename NT, class Point, class bfunc>
void test_v_poly_oracles(std::vector<Point> coeffs, bfunc phi, bfunc grad_phi, NT t_des, int facet_des) {
  typedef boost::mt19937    RNGType;
  typedef VPolytope<Point> Vpolytope;
  typedef std::tuple<NT, Point, int> result;
  Vpolytope P;
  NT tol = 1e-4;

  P = gen_cube<Vpolytope>(2, true);


  result res2 = P.curve_intersect_ipopt(0.01, 0, -1, coeffs, phi, grad_phi);
  NT t = std::get<0>(res2);

  std::cout << t << " " << t_des << std::endl;

  CHECK(std::abs(std::abs(t) - t_des) / t_des < tol);

}

template <typename NT>
void call_test_poly_oracles(char typ) {
  typedef Cartesian<NT>    Kernel;
  typedef typename Kernel::Point    Point;
  typedef std::vector<Point> pts;
  typedef std::function<NT(NT, NT, unsigned int, unsigned int)> bfunc;

  std::cout << "--- Testing intersection of 2D cube with p(t) = (t, t)" << std::endl;

  Point a0(2);
  Point a1(2);
  a1.set_coord(0, 1);
  a1.set_coord(1, 1);

  pts line_coeffs{a0, a1};

  bfunc poly_basis = [](NT t, NT t0, unsigned int j, unsigned int order) {
    return pow(t - t0, (NT) j);
  };

  bfunc poly_basis_grad = [](NT t, NT t0, unsigned int j, unsigned int order) {
    return ((NT) j) * pow(t - t0, (NT) (j - 1));
  };

  NT t_des_line = NT(1);
  int facet_des_line = 0;

  if (typ == 'H') {
    test_h_poly_oracles<NT, Point, bfunc>(line_coeffs, poly_basis, poly_basis_grad, t_des_line, facet_des_line);
  } else if (typ == 'V'){
    test_v_poly_oracles<NT, Point, bfunc>(line_coeffs, poly_basis, poly_basis_grad, t_des_line, facet_des_line);
  }

  std::cout << "--- Testing intersection of 2D cube with p(t) = (t, 2 * t^2)" << std::endl;

  Point b0(2);
  Point b1(2);
  Point b2(2);
  b1.set_coord(0, 1);
  b2.set_coord(1, 2);

  NT t_des_parabola = NT(1 / sqrt(2));
  int facet_des_parabola = 1;
  pts parabola_coeffs{b0, b1, b2};

  if (typ == 'H') {
    test_h_poly_oracles<NT, Point, bfunc>(parabola_coeffs, poly_basis, poly_basis_grad, t_des_parabola, facet_des_parabola);
  } else if (typ == 'V') {
    test_v_poly_oracles<NT, Point, bfunc>(parabola_coeffs, poly_basis, poly_basis_grad, t_des_parabola, facet_des_parabola);
  }

}

template <typename NT>
void call_benchmark_oracles() {
  typedef Cartesian<NT>    Kernel;
  typedef typename Kernel::Point    Point;
  typedef boost::mt19937    RNGType;
  typedef HPolytope<Point> Hpolytope;
  typedef std::tuple<NT, Point, int> result;
  typedef std::function<NT(NT, NT, unsigned int, unsigned int)> bfunc;
  Hpolytope P;
  NT tol = 1e-4;
  std::pair<int, int>dims = std::make_pair(2, 10);
  std::pair<int, int>orders = std::make_pair(2, 3);
  result res;
  long newton_correct = 0L, ipopt_correct = 0L, total = 0L;

  long newton_runtime = 0L;
  long ipopt_runtime = 0L;

  bfunc poly_basis = [](NT t, NT t0, unsigned int j, unsigned int order) {
    return pow(t - t0, (NT) j);
  };

  bfunc poly_basis_grad = [](NT t, NT t0, unsigned int j, unsigned int order) {
    return ((NT) j) * pow(t - t0, (NT) (j - 1));
  };


  std::vector<Point> coeffs;


  for (int dim = dims.first; dim <= dims.second; dim++) {
    Point p;
    P = gen_cube<Hpolytope>(dim, false);

    for (int order = orders.first; order <= orders.second; order++) {
      p = Point(dim);
      p.set_coord(order, order);
      coeffs.push_back(p);

      auto start = std::chrono::high_resolution_clock::now();
      res = P.curve_intersect_newton_raphson(0.01, 0, -1, coeffs, poly_basis, poly_basis_grad);
      auto stop = std::chrono::high_resolution_clock::now();

      newton_runtime += (long) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

      start = std::chrono::high_resolution_clock::now();
      P.curve_intersect_ipopt(0, 0, -1, coeffs, poly_basis, poly_basis_grad);
      stop = std::chrono::high_resolution_clock::now();
      ipopt_runtime += (long) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

      std::cout << std::endl;

      total++;
    }

  }

  std::cout << "Newton-Raphson: " << newton_runtime << " us" << std::endl;
  std::cout << "Interior-points: " << ipopt_runtime << " us" << std::endl;
}


TEST_CASE("h_poly_oracles") {
  call_test_poly_oracles<double>('H');
}

// TEST_CASE("benchmark_oracles") {
//   call_benchmark_oracles<double>();
// }

// TEST_CASE("v_poly_oracles") {
//   call_test_poly_oracles<double>('V');
// }

#endif
