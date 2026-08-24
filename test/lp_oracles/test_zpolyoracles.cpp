// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "lp_oracles/zpolyoracles.hpp"
#include <vector>
#include <tuple>

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

MT make_hexagon() {
    MT Z(3, 2);

    Z << 1, 0,
         0, 1,
         1, 1;

    return Z;
}

void test_membership() {
    MT Z = make_hexagon();
    
    std::vector<NT> row(3);
    std::vector<int> colno(3);

    CHECK(memLP_Zonotope(Z, Point(2, {0.0, 0.0})).value);
    CHECK(memLP_Zonotope(Z, Point(2, {1.0, 0.5})).value);
    CHECK(!memLP_Zonotope(Z, Point(2, {3.0, 0.0})).value);
    CHECK(!memLP_Zonotope(Z, Point(2, {0.0, 3.0})).value);
}

void test_line_intersection_vpoly() {
    MT Z = make_hexagon();

    std::vector<NT> row(4);
    std::vector<int> colno(4);

    Point p(2, {1.0, 0.0});
    Point v(2, {1.0, 0.0});

    // Single point tests.
    auto res = intersect_line_zono<NT>(Z, p, v);
    
    CHECK(res.solved);
    CHECK(res.value.first == doctest::Approx(1.0));
    CHECK(res.value.second == doctest::Approx(-3.0));
}

void test_line_intersection_zpoly() {
    MT Z(3, 2);
    
    Z << 1, 0,
         0, 1,
         1, 1;

    std::vector<NT> conv_comb(3);
    std::vector<NT> row(4);
    std::vector<int> colno(4);

    Point p(2, {0.0, 0.0});
    Point v(2, {1.0, 0.0});

    // Single point tests.
    auto res1 = intersect_line_Vpoly<NT>(Z, p, v, conv_comb.data(), true, true);
    
    auto res2 = intersect_line_Vpoly<NT>(Z, p, v, conv_comb.data(), false, true);

    CHECK(res1.solved);
    CHECK(res2.solved);
    CHECK(res1.value == doctest::Approx(-2.0));
    CHECK(res2.value == doctest::Approx(2.0));
}

TEST_CASE("test_membership") {
    test_membership();
}

TEST_CASE("line_intersection") {
    test_line_intersection_vpoly();
    test_line_intersection_zpoly();
}
