// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "lp_oracles/solve_lp.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

void make_cube(unsigned d, NT r, MT& A, VT& b) {
    A = MT::Zero(2*d, d);
    b = VT::Constant(2*d, r);

    for (unsigned i = 0; i < d; ++i) {
        A(i, i) = NT(1);
        A(d+i, i) = NT(-1);
    }
}

void test_chebychev_cube(unsigned d, NT r) {
    MT A;
    VT b;
    make_cube(d, r, A, b);

    auto res = compute_chebychev_ball<NT, Point>(A, b);
    Point px = res.value.first;
    NT rx = res.value.second;

    CHECK(res.solved);
    CHECK(rx == doctest::Approx(r));

    for (unsigned i = 0; i < d; ++i)
        CHECK(px[i] == doctest::Approx(0.0));
}

void test_chebychev_empty_cube(unsigned d) {
    MT A;
    VT b;
    make_cube(d, NT(1), A, b);

    b(0) = NT(0);
    b(d) = NT(0);

    auto res = compute_chebychev_ball<NT, Point>(A, b);
    NT rx = res.value.second;

    CHECK(res.solved);
    CHECK(rx == doctest::Approx(0.0));
}

void test_chebychev_infeasible_cube(unsigned d) {
    MT A;
    VT b;
    make_cube(d, NT(1), A, b);

    b(0) = NT(-1);
    b(d) = NT(-1);

    auto res = compute_chebychev_ball<NT, Point>(A, b);

    CHECK(!res.solved);
}

void test_identical_intersection() {
    MT V1(5, 2);

    V1 << 0, 0,
          2, 0,
          2, 2,
          0, 2,
          1, 3;

    Point direction(10);
    for (unsigned i = 0; i < 10; ++i)
        direction.set_coord(i, 0.0);

    direction.set_coord(0, NT(1));

    auto res = point_in_intersection<VT>(V1, V1, direction);
    Point p = res.value.first;
    bool empty = res.value.second;

    CHECK(res.solved);
    CHECK(!empty);
    CHECK(p[0] == doctest::Approx(0.0));
    CHECK(p[1] == doctest::Approx(0.0));
}

void test_corner_intersection() {
    MT V1(5, 2);
    MT V2(5, 2);

    V1 << 0, 0,
          2, 0,
          2, 2,
          0, 2,
          1, 3;

    V2 << 2, 2,
          4, 2,
          4, 4,
          2, 4,
          3, 5;

    Point direction(10);
    for (unsigned i = 0; i < 10; ++i)
        direction.set_coord(i, 0.0);

    direction.set_coord(2, NT(1));

    auto res = point_in_intersection<VT>(V1, V2, direction);
    Point p = res.value.first;
    bool empty = res.value.second;

    CHECK(res.solved);
    CHECK(!empty);
    CHECK(p[0] == doctest::Approx(2.0));
    CHECK(p[1] == doctest::Approx(2.0));
}

void test_empty_intersection() {
    MT V1(5, 2);
    MT V2(5, 2);

    V1 << 0, 0,
          2, 0,
          2, 2,
          0, 2,
          1, 3;

    V2 << 4, 0,
          6, 0,
          6, 2,
          4, 2,
          5, 3;

    Point direction(10);
    for (unsigned i = 0; i < 10; ++i)
        direction.set_coord(i, NT(0.0));

    auto res = point_in_intersection<VT>(V1, V2, direction);
    bool empty = res.value.second;

    CHECK(res.solved);
    CHECK(empty);
}

TEST_CASE("test_chebychev_cube") {
    test_chebychev_cube(5, NT(1));
    test_chebychev_cube(10, NT(1));
    test_chebychev_cube(15, NT(1));
    test_chebychev_empty_cube(5);
    test_chebychev_empty_cube(10);
    test_chebychev_empty_cube(15);
    test_chebychev_infeasible_cube(5);
    test_chebychev_infeasible_cube(10);
    test_chebychev_infeasible_cube(15);
}

TEST_CASE("test_intersection") {
    test_corner_intersection();
    test_identical_intersection();
    test_empty_intersection();
}
