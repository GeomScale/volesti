// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <unistd.h>
#include <string>
#include <typeinfo>

#include "doctest.h"
#include "Eigen/Eigen"

#include "ode_solvers.hpp"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/known_polytope_generators.h"


template <typename NT>
void test_euler(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;

    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    EulerODESolver<Point, NT, Hpolytope> euler_solver =
      EulerODESolver<Point, NT, Hpolytope>(0, 0.01, q, Fs, bounds{NULL});
    euler_solver.steps(1000);
    NT err=0.001;
    NT error = euler_solver.xs[0].dot(euler_solver.xs[0]);
    CHECK(error < err);
}

template <typename NT>
void test_bs(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;

    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    RichardsonExtrapolationODESolver<Point, NT, Hpolytope> bs_solver =
      RichardsonExtrapolationODESolver<Point, NT, Hpolytope>(0, 0.1, q, Fs, bounds{NULL});
    bs_solver.steps(1000);

    NT err=0.001;
    NT error = bs_solver.xs[0].dot(bs_solver.xs[0]);
    CHECK(error < err);
}


template <typename NT>
void test_rk4(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 1.0);
    pts q;
    q.push_back(q0);
    RKODESolver<Point, NT, Hpolytope> rk_solver =
      RKODESolver<Point, NT, Hpolytope>(0, 0.1, q, Fs, bounds{NULL});
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
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
    funcs Fs;
    func F = [](pts x, NT t) { return (-2.0) * x[0]; };
    Fs.push_back(F);
    Fs.push_back(F);

    // Solve in P x R for
    Hpolytope P = gen_cube<Hpolytope>(1, true);
    bounds Ks{&P, NULL};

    Point x0 = Point(1);
    Point v0 = Point(1);
    x0.set_coord(0, 0);
    v0.set_coord(0, 2.0);
    pts q{x0, v0};
    LeapfrogODESolver<Point, NT, Hpolytope> leapfrog_solver =
      LeapfrogODESolver<Point, NT, Hpolytope>(0, 0.1, q, Fs, Ks);

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
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;

    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Fs.push_back(F);
    Point x0 = Point(1);
    Point v0 = Point(1);
    x0.set_coord(0, 0);
    v0.set_coord(0, 1.0);
    pts q{x0, v0};
    LeapfrogODESolver<Point, NT, Hpolytope> leapfrog_solver =
      LeapfrogODESolver<Point, NT, Hpolytope>(0, 0.01, q, Fs, bounds{NULL, NULL});

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
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
    funcs Fs;
    bounds Ks;
    func F = [](pts xs, NT t) { return xs[0]; };
    Fs.push_back(F);

    Hpolytope P = gen_cube<Hpolytope>(1, true);
    Ks.push_back(&P);

    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    EulerODESolver<Point, NT, Hpolytope> euler_solver =
      EulerODESolver<Point, NT, Hpolytope>(0, 0.001, q, Fs, Ks);
    euler_solver.steps(1000);

    NT err=0.01;
    NT target = 1.0;
    NT error = std::abs((euler_solver.xs[0][0] - target) / target);
    CHECK(error < err);
}

template <typename NT>
void test_bs_constrained(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
    funcs Fs;
    bounds Ks;
    func F = [](pts xs, NT t) { return xs[0]; };
    Fs.push_back(F);

    Hpolytope P = gen_cube<Hpolytope>(1, true);
    Ks.push_back(&P);

    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    RichardsonExtrapolationODESolver<Point, NT, Hpolytope> bs_solver = RichardsonExtrapolationODESolver<Point, NT, Hpolytope>(0, 0.01, q, Fs, Ks);

    bs_solver.steps(1000);

    NT err=0.01;
    NT target = 1.0;
    NT error = std::abs((bs_solver.xs[0][0] - target) / target);
    CHECK(error < err);
}


template <typename NT>
void test_rk4_constrained(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
    funcs Fs;
    func F = [](pts x, NT t) { return  x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 0.2);
    pts q;
    q.push_back(q0);

    Hpolytope P = gen_cube<Hpolytope>(1, true);

    bounds Ks{&P};
    RKODESolver<Point, NT, Hpolytope> rk_solver = RKODESolver<Point, NT, Hpolytope>(0, 0.01, q, Fs, Ks);

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
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
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

    Hpolytope P = gen_cube<Hpolytope>(2, true);
    Ks.push_back(&P);

    Point q0 = Point(2);
    q0.set_coord(0, 0.2);
    q0.set_coord(1, 0.2);
    pts q;
    q.push_back(q0);
    EulerODESolver<Point, NT, Hpolytope> euler_solver = EulerODESolver<Point, NT, Hpolytope>(0, 0.1, q, Fs, Ks);
    euler_solver.steps(1000);
    CHECK(euler_solver.xs[0].dot(euler_solver.xs[0]) < 1.1);

}



template <typename NT>
void call_test_first_order() {

  std::cout << "--- Testing solution to dx / dt = -x" << std::endl;
  test_euler<NT>();
  test_rk4<NT>();
  test_bs<NT>();

  std::cout << "--- Testing solution to dx / dt = x in [-1, 1]" << std::endl;
  test_euler_constrained<NT>();
  test_rk4_constrained<NT>();
  test_bs_constrained<NT>();

  std::cout << "--- Testing solution to dx / dt = v, dv / dt = -x in [-1, 1]^2" << std::endl;
  test_euler_2d_constrained<NT>();

}

template <typename NT>
void call_test_second_order() {
  std::cout << "--- Testing solution to d^2x / dt^2 = -x" << std::endl;
  test_leapfrog<NT>();
  //
  std::cout << "--- Testing solution to d^2x / dt^2 = x in [-1, 1]" << std::endl;
  test_leapfrog_constrained<NT>();

  std::cout << "--- Testing solution to dx / dt = v, dv / dt = -x in [-1, 1]^2" << std::endl;
  test_euler_2d_constrained<NT>();

}

TEST_CASE("first_order") {
  call_test_first_order<double>();
}

TEST_CASE("second_order") {
  call_test_second_order<double>();
}

#ifndef DISABLE_NLP_ORACLES

template <typename NT>
void test_collocation(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::function<NT(NT, NT, unsigned int, unsigned int)> bfunc;
    typedef std::vector<NT> coeffs;
    typedef std::vector<func> funcs;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;

    funcs Fs;
    func F = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(F);
    Point q0 = Point(1);
    q0.set_coord(0, 1.0);
    pts q;
    q.push_back(q0);

    bfunc phi = [](NT t, NT t0, unsigned int j, unsigned int order) {
      return pow(t - t0, (NT) j);
    };

    bfunc grad_phi = [](NT t, NT t0, unsigned int j, unsigned int order) {
      return ((NT) j) * pow(t - t0, (NT) (j - 1));
    };

    // Trapezoidal collocation
    coeffs cs{0.0, 0.0, 1.0};


    CollocationODESolver<Point, NT, Hpolytope, bfunc> c_solver =
      CollocationODESolver<Point, NT, Hpolytope, bfunc>
      (0, 1.0, q, Fs, bounds{NULL}, cs, phi, grad_phi);
    c_solver.steps(100);
    NT err=0.001;
    NT error = c_solver.xs[0].dot(c_solver.xs[0]);
    CHECK(error < err);
}

template <typename NT>
void test_collocation_constrained(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef std::vector<Point> pts;
    typedef std::function<Point(pts, NT)> func;
    typedef std::vector<func> funcs;
    typedef std::function<NT(NT, NT, unsigned int, unsigned int)> bfunc;
    typedef std::vector<NT> coeffs;
    typedef HPolytope<Point>  Hpolytope;
    typedef std::vector<Hpolytope*> bounds;
    funcs Fs;
    bounds Ks;
    func F = [](pts xs, NT t) { return xs[0]; };
    Fs.push_back(F);

    Hpolytope P = gen_cube<Hpolytope>(1, false);
    Ks.push_back(&P);

    bfunc phi = [](NT t, NT t0, unsigned int j, unsigned int order) {
      return pow(t - t0, (NT) j);
    };

    bfunc grad_phi = [](NT t, NT t0, unsigned int j, unsigned int order) {
      return ((NT) j) * pow(t - t0, (NT) (j - 1));
    };

    // Trapezoidal collocation
    coeffs cs{0.0, 0.0, 1.0};

    Point q0 = Point(1);
    q0.set_coord(0, 0.5);
    pts q;
    q.push_back(q0);
    CollocationODESolver<Point, NT, Hpolytope, bfunc> c_solver =
      CollocationODESolver<Point, NT, Hpolytope, bfunc>
      (0, 0.05, q, Fs, Ks, cs, phi, grad_phi);
    c_solver.steps(1000);


    NT err=0.1;
    NT target = 1.0;
    NT error = std::abs((c_solver.xs[0].dot(c_solver.xs[0]) - target) / target);
    CHECK(error < err);
}

template <typename NT>
void call_test_collocation() {

  std::cout << "--- Testing solution to dx / dt = -x w/ collocation" << std::endl;
  test_collocation<NT>();

  std::cout << "--- Testing solution to dx / dt = x in [-1, 1] w/ collocation" << std::endl;
  test_collocation_constrained<NT>();

}

TEST_CASE("collocation") {
  call_test_collocation<double>();
}

#endif
