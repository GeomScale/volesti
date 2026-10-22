// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

// Edge cases of H-polytope preprocessing and of the sequence-of-balls volume
// algorithm. Built twice: with lp_solve and with DISABLE_LPSOLVE, since the
// inner-ball failure path differs between the two.

#include "doctest.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef HPolytope<Point> Hpolytope;
typedef HPolytope<Point, Eigen::SparseMatrix<NT, Eigen::RowMajor>> SparseHpolytope;

// The box [lo, hi]^2 as {x : A x <= b}, optionally with an extra zero row 0 <= b0
template <typename Polytope>
Polytope box2(NT lo, NT hi, bool zero_row = false, NT b0 = NT(0))
{
    int m = zero_row ? 5 : 4;
    MT A = MT::Zero(m, 2);
    VT b(m);
    A << MT::Identity(2, 2), -MT::Identity(2, 2), MT::Zero(m - 4, 2);
    b.head(4) << hi, hi, -lo, -lo;
    if (zero_row) b(4) = b0;
    return Polytope(2, A, b);
}

template <typename Polytope>
void check_ball(Polytope P, NT radius, NT c0, NT c1)
{
    auto ball = P.ComputeInnerBall();
    CHECK(ball.second == doctest::Approx(radius).epsilon(1e-6));
    CHECK(ball.first[0] == doctest::Approx(c0).epsilon(1e-6));
    CHECK(ball.first[1] == doctest::Approx(c1).epsilon(1e-6));
}

template <typename Polytope>
void call_test_zero_rows()
{
    // A zero row 0 <= b_i is not a facet: it must not cap the radius at b_i
    check_ball(box2<Polytope>(-1, 1), 1, 0, 0);
    check_ball(box2<Polytope>(-1, 1, true, 0), 1, 0, 0);
    check_ball(box2<Polytope>(-1, 1, true, 0.5), 1, 0, 0);
    // ... and 0 <= b_i with b_i < 0 makes the polytope empty
    auto ball = box2<Polytope>(-1, 1, true, -1).ComputeInnerBall();
    CHECK(ball.second < 0);
}

template <typename Polytope>
void call_test_inner_ball_failure()
{
    // Empty polytope: the radius must be negative and the center must have
    // the polytope's dimension (callers check second < 0)
    MT A(4, 2);
    A << 1, 0, 0, 1, -1, 0, 0, -1;
    VT b(4);
    b << 1, 1, -2, 1;   // x <= 1 and x >= 2
    Polytope P(2, A, b);
    auto ball = P.ComputeInnerBall();
    CHECK(ball.second < 0);
    CHECK(ball.first.dimension() == 2);
    // The failure is cached, as on the lp_solve path
    CHECK(P.ComputeInnerBall().second < 0);
    // Flat polytope (x = 1): no inscribed ball of positive radius
    b << 1, 1, -1, 1;
    Polytope Q(2, A, b);
    CHECK(!(Q.ComputeInnerBall().second > 0));
}

template <typename Polytope>
void call_test_set_interior_point()
{
    Polytope P = box2<Polytope>(-1, 1, true, 0);
    Point p(2);
    p.set_coord(0, 0.5);
    P.set_interior_point(p);
    // the zero row does not count: the distance to the nearest facet is 0.5
    CHECK(P.InnerBall().second == doctest::Approx(0.5));

    CHECK_THROWS_AS(P.set_interior_point(Point(3)), std::invalid_argument);
    Point q(2);
    q.set_coord(0, std::numeric_limits<NT>::quiet_NaN());
    CHECK_THROWS_AS(P.set_interior_point(q), std::invalid_argument);
    q.set_coord(0, 2);   // outside
    CHECK_THROWS_AS(P.set_interior_point(q), std::invalid_argument);
    q.set_coord(0, 1);   // on the boundary
    CHECK_THROWS_AS(P.set_interior_point(q), std::invalid_argument);
}

void call_test_sob_small_radius()
{
    // Regression: for an inner radius < 1 the first ball index was truncated
    // towards zero, and when it met the last index the "volume" was that of
    // the inner ball alone (pi * 0.5005^2 = 0.787 here, area 1.002001)
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    Hpolytope P = box2<Hpolytope>(-0.001, 1);
    NT exact = 1.001 * 1.001;
    NT volume = volume_sequence_of_balls<CDHRWalk, RNGType>(P, 0.1, 1);
    std::cout << "Computed volume " << volume << ", exact " << exact << std::endl;
    CHECK(std::abs(volume - exact) / exact < 0.1);
}

void call_test_sob_inner_ball_failure()
{
    // An empty polytope has no inner ball: SoB throws instead of returning -1
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    MT A(4, 2);
    A << 1, 0, 0, 1, -1, 0, 0, -1;
    VT b(4);
    b << 1, 1, -2, 1;   // x <= 1 and x >= 2
    Hpolytope P(2, A, b);
    CHECK_THROWS_AS((volume_sequence_of_balls<CDHRWalk, RNGType>(P, 0.1, 1)),
                    std::runtime_error);
}

TEST_CASE("zero_rows") {
    call_test_zero_rows<Hpolytope>();
    call_test_zero_rows<SparseHpolytope>();
}

TEST_CASE("inner_ball_failure") {
    call_test_inner_ball_failure<Hpolytope>();
    call_test_inner_ball_failure<SparseHpolytope>();
}

TEST_CASE("set_interior_point") {
    call_test_set_interior_point<Hpolytope>();
    call_test_set_interior_point<SparseHpolytope>();
}

TEST_CASE("sob_small_radius") {
    call_test_sob_small_radius();
}

TEST_CASE("sob_inner_ball_failure") {
    call_test_sob_inner_ball_failure();
}
