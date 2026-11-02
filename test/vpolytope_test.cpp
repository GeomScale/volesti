// VolEsti (volume computation and sampling library)

// Licensed under GNU LGPL.3, see LICENCE file

// Regression tests for the special member functions (copy/move construction
// and assignment) of VPolytope. VPolytope owns scratch buffers that are passed
// to the lp_solve based oracles, so a copy or a move that leaves them dangling
// or null breaks every later oracle call. Run under AddressSanitizer
// (-fsanitize=address) to also detect use-after-free and leaks.

#include "doctest.h"

#include <utility>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/vpolytope.h"
#include "generators/known_polytope_generators.h"

typedef Cartesian<double>             Kernel;
typedef Kernel::Point                 Point;
typedef VPolytope<Point>              Vpolytope;

// The cross polytope conv{+-e_i} contains the origin, so the ray from the
// origin along e_0 leaves it at distance 1 (and at -1 in the opposite direction).
static std::pair<double, double> shoot_ray(Vpolytope const& P)
{
    unsigned int d = P.dimension();
    Point origin(d), v(d);
    for (unsigned int i = 0; i < d; ++i) {
        origin.set_coord(i, 0.0);
        v.set_coord(i, i == 0 ? 1.0 : 0.0);
    }
    return P.line_intersect(origin, v);
}

static void check_oracles(Vpolytope const& P)
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

TEST_CASE("vpolytope_copy_constructor") {
    Vpolytope P = generate_cross<Vpolytope>(3, true);
    Vpolytope Q(P);
    check_oracles(Q);
    check_oracles(P);
}

TEST_CASE("vpolytope_copy_assignment") {
    Vpolytope P = generate_cross<Vpolytope>(3, true);
    Vpolytope Q = generate_cross<Vpolytope>(3, true);
    Q = P;
    check_oracles(Q);
    check_oracles(P);
}

TEST_CASE("vpolytope_move_constructor") {
    Vpolytope P = generate_cross<Vpolytope>(3, true);
    Vpolytope Q(std::move(P));
    check_oracles(Q);
}

TEST_CASE("vpolytope_move_assignment") {
    Vpolytope Q = generate_cross<Vpolytope>(3, true);
    Q = generate_cross<Vpolytope>(3, true);
    check_oracles(Q);
}

TEST_CASE("vpolytope_default_then_move_assignment") {
    // the pattern used throughout the volume tests
    Vpolytope P;
    P = generate_cross<Vpolytope>(3, true);
    check_oracles(P);
}

TEST_CASE("vpolytope_copy_and_move_preserve_inner_ball") {
    Vpolytope P = generate_cross<Vpolytope>(3, true);
    P.set_InnerBall(std::make_pair(Point(3), 0.5));

    Vpolytope copy_constructed(P);
    CHECK(copy_constructed.InnerBall().second == doctest::Approx(0.5));

    Vpolytope copy_assigned;
    copy_assigned = P;
    CHECK(copy_assigned.InnerBall().second == doctest::Approx(0.5));

    Vpolytope move_constructed(std::move(P));
    CHECK(move_constructed.InnerBall().second == doctest::Approx(0.5));

    Vpolytope move_assigned;
    move_assigned = std::move(move_constructed);
    CHECK(move_assigned.InnerBall().second == doctest::Approx(0.5));
}
