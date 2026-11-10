// VolEsti (volume computation and sampling library)

// Licensed under GNU LGPL.3, see LICENCE file

// Regression tests for the special member functions (copy/move construction
// and assignment) of Zonotope. Zonotope owns scratch buffers that are passed
// to the lp_solve based oracles, so a copy or a move that leaves them dangling
// or null breaks every later oracle call. Run under AddressSanitizer
// (-fsanitize=address) to also detect use-after-free and leaks.

#include "doctest.h"

#include <utility>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/zpolytope.h"

typedef Cartesian<double>             Kernel;
typedef Kernel::Point                 Point;
typedef Zonotope<Point>               zonotope;
typedef zonotope::MT                  MT;
typedef zonotope::VT                  VT;

// Helper to construct a hypercube zonotope [-1, 1]^d from standard basis vectors
static zonotope make_cube_zonotope(unsigned int d)
{
    MT V = MT::Identity(d, d);
    VT b = VT::Ones(d);
    return zonotope(d, V, b);
}

static std::pair<double, double> shoot_ray(zonotope const& P)
{
    unsigned int d = P.dimension();
    Point origin(d), v(d);
    for (unsigned int i = 0; i < d; ++i) {
        origin.set_coord(i, 0.0);
        v.set_coord(i, i == 0 ? 1.0 : 0.0);
    }
    return P.line_intersect(origin, v);
}

static void check_oracles(zonotope const& P)
{
    std::pair<double, double> res = shoot_ray(P);
    CHECK(res.first == doctest::Approx(1.0));
    CHECK(res.second == doctest::Approx(-1.0));

    unsigned int d = P.dimension();
    Point inside(d), outside(d);
    for (unsigned int i = 0; i < d; ++i) {
        inside.set_coord(i, 0.0);
        outside.set_coord(i, 2.0);
    }
    CHECK(P.is_in(inside) == -1);
    CHECK(P.is_in(outside) == 0);
}

TEST_CASE("zonotope_copy_constructor") {
    zonotope P = make_cube_zonotope(3);
    zonotope Q(P);
    check_oracles(Q);
    check_oracles(P);
}

TEST_CASE("zonotope_copy_assignment") {
    zonotope P = make_cube_zonotope(3);
    zonotope Q = make_cube_zonotope(3);
    Q = P;
    check_oracles(Q);
    check_oracles(P);
}

TEST_CASE("zonotope_move_constructor") {
    zonotope P = make_cube_zonotope(3);
    zonotope Q(std::move(P));
    check_oracles(Q);
}

TEST_CASE("zonotope_move_assignment") {
    zonotope Q = make_cube_zonotope(3);
    Q = make_cube_zonotope(3);
    check_oracles(Q);
}

TEST_CASE("zonotope_default_then_move_assignment") {
    // the pattern used throughout the volume tests
    zonotope P;
    P = make_cube_zonotope(3);
    check_oracles(P);
}

TEST_CASE("zonotope_copy_and_move_preserve_inner_ball") {
    zonotope P = make_cube_zonotope(3);
    P.set_InnerBall(std::make_pair(Point(3), 0.5));

    zonotope copy_constructed(P);
    CHECK(copy_constructed.InnerBall().second == doctest::Approx(0.5));

    zonotope copy_assigned;
    copy_assigned = P;
    CHECK(copy_assigned.InnerBall().second == doctest::Approx(0.5));

    zonotope move_constructed(std::move(P));
    CHECK(move_constructed.InnerBall().second == doctest::Approx(0.5));

    zonotope move_assigned;
    move_assigned = std::move(move_constructed);
    CHECK(move_assigned.InnerBall().second == doctest::Approx(0.5));
}
