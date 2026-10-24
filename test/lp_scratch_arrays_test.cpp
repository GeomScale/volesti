// VolEsti (volume computation and sampling library)

// Licensed under GNU LGPL.3, see LICENCE file

// Tests for the lp_solve based helpers in lp_oracles/solve_lp.h. Besides their
// results, these tests are meant to be run under LeakSanitizer
// (-fsanitize=address): the helpers used to leak their scratch arrays on every call.

#include "doctest.h"

#include <Eigen/Eigen>

#include "cartesian_geom/cartesian_kernel.h"
#include "lp_oracles/solve_lp.h"

typedef double                                              NT;
typedef Cartesian<NT>                                       Kernel;
typedef Kernel::Point                                       Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;

// the vertices of the square [lo, lo + side]^2, one vertex per row
static MT square_vertices(NT lo, NT side)
{
    MT V(4, 2);
    V << lo,        lo,
         lo + side, lo,
         lo,        lo + side,
         lo + side, lo + side;
    return V;
}

TEST_CASE("chebychev_ball_of_cube") {
    // the cube [-1, 1]^3 written as A x <= b
    unsigned int const d = 3;
    MT A(2 * d, d);
    A.setZero();
    VT b = VT::Ones(2 * d);
    for (unsigned int i = 0; i < d; ++i) {
        A(i, i) = 1.0;
        A(d + i, i) = -1.0;
    }

    // called repeatedly, as ComputeInnerBall() is during volume computation
    for (int rep = 0; rep < 5; ++rep) {
        std::pair<Point, NT> ball = ComputeChebychevBall<NT, Point>(A, b);
        CHECK(ball.second == doctest::Approx(1.0));
        for (unsigned int i = 0; i < d; ++i) {
            CHECK(ball.first[i] == doctest::Approx(0.0).epsilon(1e-9));
        }
    }
}

TEST_CASE("point_in_intersection_of_squares") {
    MT V1 = square_vertices(0.0, 2.0);   // [0, 2]^2
    MT V2 = square_vertices(1.0, 2.0);   // [1, 3]^2, overlaps V1 in [1, 2]^2
    Point direction(8);                  // one objective coefficient per vertex of both squares
    for (int j = 0; j < 8; ++j) direction.set_coord(j, 1.0);

    for (int rep = 0; rep < 5; ++rep) {
        bool empty = true;
        Point p = PointInIntersection<VT>(V1, V2, direction, empty);
        CHECK(!empty);
        for (int i = 0; i < 2; ++i) {
            CHECK(p[i] >= 1.0 - 1e-9);
            CHECK(p[i] <= 2.0 + 1e-9);
        }
    }
}

TEST_CASE("point_in_intersection_of_disjoint_squares") {
    MT V1 = square_vertices(0.0, 1.0);   // [0, 1]^2
    MT V2 = square_vertices(5.0, 1.0);   // [5, 6]^2
    Point direction(8);
    for (int j = 0; j < 8; ++j) direction.set_coord(j, 1.0);

    for (int rep = 0; rep < 5; ++rep) {
        bool empty = false;
        PointInIntersection<VT>(V1, V2, direction, empty);
        CHECK(empty);
    }
}
