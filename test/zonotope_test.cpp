// VolEsti (volume computation and sampling library)

// Licensed under GNU LGPL.3, see LICENCE file

// Regression tests for the special member functions (copy/move construction
// and assignment) of Zonotope. Zonotope owns scratch buffers that are passed to
// the lp_solve based oracles, so a copy or a move that leaves them dangling or
// null breaks every later oracle call. Run under AddressSanitizer
// (-fsanitize=address) to also detect use-after-free and leaks.

#include "doctest.h"

#include <chrono>
#include <cmath>
#include <utility>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/zpolytope.h"
#include "generators/z_polytopes_generators.h"

typedef Cartesian<double>             Kernel;
typedef Kernel::Point                 Point;
typedef Zonotope<Point>               Zono;
typedef boost::mt19937                RNGType;

static const int dim = 3;
static const int num_generators = 6;

// a deterministic zonotope; it is centrally symmetric around the origin
static Zono make_zonotope()
{
    return gen_zonotope_uniform<Zono, RNGType>(dim, num_generators, 7);
}

static std::pair<double, double> shoot_ray(Zono const& Z)
{
    Point origin(dim), v(dim);
    for (int i = 0; i < dim; ++i) {
        origin.set_coord(i, 0.0);
        v.set_coord(i, i == 0 ? 1.0 : 0.0);
    }
    return Z.line_intersect(origin, v);
}

static void check_oracles(Zono const& Z)
{
    // reference values from a freshly constructed, never copied or moved zonotope
    Zono const ref = make_zonotope();
    std::pair<double, double> expected = shoot_ray(ref);
    std::pair<double, double> res = shoot_ray(Z);
    CHECK(res.first > 0.0);
    CHECK(res.first == doctest::Approx(expected.first));
    CHECK(res.second == doctest::Approx(expected.second));

    Point inside(dim), outside(dim);
    for (int i = 0; i < dim; ++i) {
        inside.set_coord(i, 0.0);
        outside.set_coord(i, 1e6);
    }
    CHECK(Z.is_in(inside) == -1);
    CHECK(Z.is_in(outside) == 0);
}

TEST_CASE("zonotope_copy_constructor") {
    Zono Z = make_zonotope();
    Zono W(Z);
    check_oracles(W);
    check_oracles(Z);
}

TEST_CASE("zonotope_copy_assignment") {
    Zono Z = make_zonotope();
    Zono W = make_zonotope();
    W = Z;
    check_oracles(W);
    check_oracles(Z);
}

TEST_CASE("zonotope_move_constructor") {
    Zono Z = make_zonotope();
    Zono W(std::move(Z));
    check_oracles(W);
}

TEST_CASE("zonotope_move_assignment") {
    Zono W = make_zonotope();
    W = make_zonotope();
    check_oracles(W);
}

TEST_CASE("zonotope_default_then_move_assignment") {
    Zono Z;
    Z = make_zonotope();
    check_oracles(Z);
}

TEST_CASE("zonotope_copy_and_move_preserve_data") {
    Zono const ref = make_zonotope();
    Zono Z = make_zonotope();
    Z.set_InnerBall(std::make_pair(Point(dim), 0.5));

    Zono copy_constructed(Z);
    Zono copy_assigned;
    copy_assigned = Z;
    Zono move_constructed(std::move(Z));
    Zono move_assigned;
    move_assigned = std::move(move_constructed);

    Zono const* results[] = {&copy_constructed, &copy_assigned, &move_assigned};
    for (Zono const* R : results) {
        CHECK(R->InnerBall().second == doctest::Approx(0.5));
        // the data computed by compute_eigenvectors is carried over
        REQUIRE(R->get_Q0().rows() == ref.get_Q0().rows());
        REQUIRE(R->get_Q0().cols() == ref.get_Q0().cols());
        CHECK((R->get_Q0() - ref.get_Q0()).norm() < 1e-12);
        CHECK((R->get_mat() - ref.get_mat()).norm() < 1e-12);
    }
}
