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
#include "ode_solvers.h"
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include "Eigen/Eigen"
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

template <typename NT>
void test_euler(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    EulerODESolver<Point, NT, Vpolytope> euler_solver = EulerODESolver<Point, NT, Vpolytope>(0, 0.01, q, Fs);
    euler_solver.steps(1000);
    NT err=0.001;
    NT error = euler_solver.xs[0].dot(euler_solver.xs[0]);
    CHECK(error < err);
}


template <typename NT>
void test_rk4(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 1.0);
    pts q;
    q.push_back(q0);
    RKODESolver<Point, NT, Vpolytope> rk_solver = RKODESolver<Point, NT, Vpolytope>(0, 0.1, q, Fs);
    rk_solver.steps(1000);

    NT err=0.001;
    NT error = rk_solver.xs[0].dot(rk_solver.xs[0]);
    CHECK(error < err);
}


template <typename NT>
void test_leapfrog_constrained(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef std::vector<Vpolytope*> bounds;
    funcs Fs;
    func F = [](pts x, NT t) { return (-2.0) * x[0]; };
    Fs.push_back(F);
    Fs.push_back(F);

    // Solve in P x R for
    Vpolytope P = gen_cube<Vpolytope>(1, true);
    bounds Ks{&P, NULL};

    Point x0 = Point(1);
    Point v0 = Point(1);
    x0.set_coord(0, 0);
    v0.set_coord(0, 2.0);
    pts q{x0, v0};
    LeapfrogODESolver<Point, NT, Vpolytope> leapfrog_solver = LeapfrogODESolver<Point, NT, Vpolytope>(0, 0.01, q, Fs, Ks);

    for (int i = 0; i < 1000; i++) {
      leapfrog_solver.step();
      CHECK(leapfrog_solver.xs[0].dot(leapfrog_solver.xs[0]) < 1.1);
    }

}

template <typename NT>
void test_leapfrog(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Fs.push_back(F);
    Point x0 = Point(1);
    Point v0 = Point(1);
    x0.set_coord(0, 0);
    v0.set_coord(0, 1.0);
    pts q{x0, v0};
    LeapfrogODESolver<Point, NT, Vpolytope> leapfrog_solver = LeapfrogODESolver<Point, NT, Vpolytope>(0, 0.01, q, Fs);

    for (int i = 0; i < 1000; i++) {
      leapfrog_solver.step();
      CHECK(leapfrog_solver.xs[0].dot(leapfrog_solver.xs[0]) < 1.1);
    }

}


template <typename NT>
void test_euler_constrained(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef std::vector<Vpolytope*> bounds;
    funcs Fs;
    bounds Ks;
    func F = [](pts xs, NT t) { return xs[0]; };
    Fs.push_back(F);

    Vpolytope P = gen_cube<Vpolytope>(1, true);
    Ks.push_back(&P);

    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    EulerODESolver<Point, NT, Vpolytope> euler_solver = EulerODESolver<Point, NT, Vpolytope>(0, 0.001, q, Fs, Ks);
    euler_solver.steps(1000);

    NT err=0.01;
    NT target = 1.0;
    NT error = std::abs((euler_solver.xs[0][0] - target) / target);
    CHECK(error < err);
}



template <typename NT>
void test_rk4_constrained(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef std::vector<Vpolytope*> bounds;
    funcs Fs;
    func F = [](pts x, NT t) { return  x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 0.2);
    pts q;
    q.push_back(q0);

    Vpolytope P = gen_cube<Vpolytope>(1, true);

    bounds Ks{&P};
    RKODESolver<Point, NT, Vpolytope> rk_solver = RKODESolver<Point, NT, Vpolytope>(0, 0.01, q, Fs, Ks);

    rk_solver.steps(1000);

    NT err=0.01;
    NT target = 1.0;
    NT error = std::abs((rk_solver.xs[0][0] - target) / target);
    CHECK(error < err);
}


template <typename NT>
void test_euler_2d_constrained(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef std::vector<Vpolytope*> bounds;
    funcs Fs;
    bounds Ks;
    func F = [](pts xs, NT t) {
        Point y(xs[0].dimension());
        // y.set_coord(0, 2 * (xs[0][0] - 1 / 3 * pow(xs[0][0], 3) - xs[0][1]));
        // y.set_coord(1, (1 / 2) * xs[0][0]);
        y.set_coord(0, (1) * xs[0][1]);
        y.set_coord(1, (-1.0) * xs[0][0]);

        return y;
     };

    Fs.push_back(F);

    Vpolytope P = gen_cube<Vpolytope>(2, true);
    Ks.push_back(&P);

    Point q0 = Point(2);
    q0.set_coord(0, 0.2);
    q0.set_coord(1, 0.2);
    pts q;
    q.push_back(q0);
    EulerODESolver<Point, NT, Vpolytope> euler_solver = EulerODESolver<Point, NT, Vpolytope>(0, 0.1, q, Fs, Ks);
    euler_solver.steps(1000);
    CHECK(euler_solver.xs[0].dot(euler_solver.xs[0]) < 1.1);

}



template <typename NT>
void call_test_euler() {

  std::cout << "--- Testing solution to dx / dt = -x" << std::endl;
  test_euler<NT>();
  test_rk4<NT>();

  std::cout << "--- Testing solution to dx / dt = x in [-1, 1]" << std::endl;
  test_euler_constrained<NT>();
  test_rk4_constrained<NT>();

  std::cout << "--- Testing solution to dx / dt = v, dv / dt = -x in [-1, 1]^2" << std::endl;
  test_euler_2d_constrained<NT>();

}


template <typename NT>
void call_test_leapfrog() {
  std::cout << "--- Testing solution to d^2x / dt^2 = -x" << std::endl;
  test_leapfrog<NT>();
  //
  std::cout << "--- Testing solution to d^2x / dt^2 = x in [-1, 1]" << std::endl;
  test_leapfrog_constrained<NT>();
  // //
  std::cout << "--- Testing solution to dx / dt = v, dv / dt = -x in [-1, 1]^2" << std::endl;
  test_euler_2d_constrained<NT>();

}

TEST_CASE("euler") {
  call_test_euler<double>();
}

TEST_CASE("leapfrog") {
  call_test_leapfrog<double>();
}
