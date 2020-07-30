// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DISABLE_NLP_ORACLES

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include <string>
#include <typeinfo>
#include <chrono>

#include "Eigen/Eigen"
#include "doctest.h"

#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/known_polytope_generators.h"
#include "nlp_oracles/nlp_oracles.hpp"

template <typename NT, class Point, class bfunc>
void test_h_poly_oracles(std::vector<Point> coeffs, bfunc phi, bfunc grad_phi, NT t_des, int facet_des) {
  typedef boost::mt19937    RNGType;
  typedef HPolytope<Point> Hpolytope;
  typedef std::tuple<NT, Point, int> result;
  Hpolytope P;
  NT tol = 1e-4;

  P = gen_cube<Hpolytope>(2, false);
  NewtonRaphsonHPolyoracle<Hpolytope, bfunc> nr_oracle;
  IpoptHPolyoracle<Hpolytope, bfunc> ipopt_oracle;
  MPSolveHPolyoracle<Hpolytope, bfunc> mpsolve_oracle;

  result res = P.curve_intersect(0.01, 0, -1, coeffs, phi, grad_phi, nr_oracle);
  NT t = std::get<0>(res);
  int facet = std::get<2>(res);

  // CHECK(facet == facet_des);
  CHECK(std::abs(t - t_des) / t_des < tol);

  result res2 = P.curve_intersect(0.0, 0, -1, coeffs, phi, grad_phi, ipopt_oracle);

  t = std::get<0>(res2);
  facet = std::get<2>(res2);
  CHECK(std::abs(t - t_des) / t_des < tol);

  res2 = P.curve_intersect(0, 0, -1, coeffs, phi, grad_phi, mpsolve_oracle);

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
  IpoptVPolyoracle<Vpolytope, bfunc> ipopt_oracle;


  result res2 = P.curve_intersect(0.01, 0, -1, coeffs, phi, grad_phi, ipopt_oracle);
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

  NewtonRaphsonHPolyoracle<Hpolytope, bfunc> nr_oracle;
  IpoptHPolyoracle<Hpolytope, bfunc> ipopt_oracle;

  for (int dim = dims.first; dim <= dims.second; dim++) {
    Point p;
    P = gen_cube<Hpolytope>(dim, false);

    for (int order = orders.first; order <= orders.second; order++) {
      p = Point(dim);
      p.set_coord(order, order);
      coeffs.push_back(p);

      auto start = std::chrono::high_resolution_clock::now();
      res = P.curve_intersect(0.01, 0, -1, coeffs, poly_basis, poly_basis_grad, nr_oracle);
      auto stop = std::chrono::high_resolution_clock::now();

      newton_runtime += (long) std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

      start = std::chrono::high_resolution_clock::now();
      P.curve_intersect(0, 0, -1, coeffs, poly_basis, poly_basis_grad, ipopt_oracle);
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
