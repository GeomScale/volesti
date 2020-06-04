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
#include "volume.h"
#include "known_polytope_generators.h"
#include <string>
#include <typeinfo>
#include "samplers.h"
#include "doctest.h"

template <typename NT, class Point, class bfunc>
void test_h_poly_oracles(std::vector<Point> coeffs, bfunc phi, bfunc grad_phi) {
  typedef boost::mt19937    RNGType;
  typedef HPolytope<Point> Hpolytope;
  typedef std::vector<std::pair<NT, int>> results;
  Hpolytope P;

  P = gen_cube<Hpolytope>(2, false);

  results res = P.curve_intersect(0.01, 0, coeffs, phi, grad_phi);

  std::cout << "t = " << res[0].first << " facet = " << res[0].second << std::endl;


}

template <typename NT>
void call_test_h_poly_oracles() {
  typedef Cartesian<NT>    Kernel;
  typedef typename Kernel::Point    Point;
  typedef std::vector<Point> pts;
  typedef std::function<NT(NT, NT, unsigned int, unsigned int)> bfunc;

  std::cout << "--- Testing intersection of 2D H-cube with p(t) = (t, t)" << std::endl;

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

  test_h_poly_oracles<NT, Point, bfunc>(line_coeffs, poly_basis, poly_basis_grad);

  std::cout << "--- Testing intersection of 2D H-cube with p(t) = (t, 2 * t^2)" << std::endl;

  Point b0(2);
  Point b1(2);
  Point b2(2);
  b1.set_coord(0, 1);
  b2.set_coord(1, 2);

  pts parabola_coeffs{b0, b1, b2};

  test_h_poly_oracles<NT, Point, bfunc>(parabola_coeffs, poly_basis, poly_basis_grad);


}


TEST_CASE("h_poly_oracles") {
  call_test_h_poly_oracles<double>();
}
