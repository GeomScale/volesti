// VolEsti (volume computation and sampling library)
// Tests for small utility files with 0% coverage
#include "doctest.h"
#include "generators/known_polytope_generators.h"
#include "generators/convex_bodies_generator.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/convex_body.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/max_inscribed_ball.hpp"

// =============================================================================
// Known polytope generators
// =============================================================================

template <typename NT> void call_test_known_generators() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point> Vpolytope;
    unsigned int dim = 2;

    SUBCASE("cube H") {
        Hpolytope P = generate_cube<Hpolytope>(dim, false);
        CHECK(P.dimension() == dim);
        P.ComputeInnerBall();
        Point origin(dim);
        CHECK(P.is_in(origin) == -1);
    }
    SUBCASE("cube custom scale") {
        Hpolytope P = generate_cube<Hpolytope>(dim, false, NT(2.0));
        CHECK(P.dimension() == dim);
    }
    SUBCASE("cube V") {
        Vpolytope P = generate_cube<Vpolytope>(dim, true);
        CHECK(P.dimension() == dim);
        Point origin(dim);
        P.ComputeInnerBall();
        CHECK(P.is_in(origin) == -1);
    }
    SUBCASE("cross H") {
        Hpolytope P = generate_cross<Hpolytope>(dim, false);
        CHECK(P.dimension() == dim);
        Point origin(dim);
        P.ComputeInnerBall();
        CHECK(P.is_in(origin) == -1);
    }
    SUBCASE("simplex H") {
        Hpolytope P = generate_simplex<Hpolytope>(dim, false);
        CHECK(P.dimension() == dim);
        Point origin(dim);
        CHECK(P.is_in(origin) == -1);
    }
    SUBCASE("prod_simplex") {
        Hpolytope P = generate_prod_simplex<Hpolytope>(dim);
        CHECK(P.dimension() == 2 * dim);
    }
    SUBCASE("skinny_cube") {
        Hpolytope P = generate_skinny_cube<Hpolytope>(dim);
        CHECK(P.dimension() == dim);
    }
    SUBCASE("birkhoff") {
        Hpolytope P = generate_birkhoff<Hpolytope>(2);
        CHECK(P.dimension() == 1);
        Hpolytope P3 = generate_birkhoff<Hpolytope>(3);
        CHECK(P3.dimension() == 4);
    }
}

TEST_CASE("known_polytope_generators") {
    call_test_known_generators<double>();
}

// =============================================================================
// ConvexBody and generator tests
// =============================================================================

template <typename NT> void call_test_convex_body() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef ConvexBody<Point> CB;
    int dim = 2;

    SUBCASE("unit ball") {
        CB B = generate_unit_ball<CB>(dim);
        CHECK(B.dimension() == dim);
        Point origin(dim);
        CHECK(B.is_in(origin) == -1);
    }
    SUBCASE("unit ball intersect hyperplane") {
        CB B = generate_unit_ball_intersect_hyperplane<CB>(dim);
        CHECK(B.dimension() == dim);
        Point origin(dim);
        CHECK(B.is_in(origin) == -1);
    }
    SUBCASE("unit ball intersect logsumexp") {
        CB B = generate_unit_ball_intersect_logsumexp<CB>(2);
        CHECK(B.dimension() == 2);
    }
}

TEST_CASE("convex_body") {
    call_test_convex_body<double>();
}

// =============================================================================
// More Poset tests
// =============================================================================

#include "misc/poset.h"

TEST_CASE("Poset transitive chain") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0,1},{1,2},{2,3}};
    Poset p(4, rel);
    CHECK(p.num_elem() == 4);
    CHECK(p.num_relations() == 3);
}
