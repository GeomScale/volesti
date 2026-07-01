// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "doctest.h"

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/convex_body.h"
#include "generators/convex_bodies_generator.h"

template <typename NT>
void test_generate_unit_ball() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef ConvexBody<Point> CB;

    unsigned int dim = 2;
    CB ball = generate_unit_ball<CB>(dim);

    CHECK(ball.dimension() == dim);

    // Origin should be inside the unit ball
    Point origin(dim);
    CHECK(ball.is_in(origin) == -1);

    // Point far outside (norm > 1) should be outside
    Point outside(dim);
    outside.set_coord(0, 2.0);
    CHECK(ball.is_in(outside) == 0);

    // Point just inside (norm < 1) should be inside
    Point inside(dim);
    inside.set_coord(0, 0.5);
    CHECK(ball.is_in(inside) == -1);

    // Point on the boundary (norm = 1): constraint g(p) = x·x - 1 = 0.
    // is_in checks g(p) > 0, so 0 > 0 is false and it returns -1 (inside).
    Point boundary(dim);
    boundary.set_coord(0, 1.0);
    CHECK(ball.is_in(boundary) == -1);
}

template <typename NT>
void test_generate_unit_ball_intersect_hyperplane() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef ConvexBody<Point> CB;

    unsigned int dim = 2;
    CB body = generate_unit_ball_intersect_hyperplane<CB>(dim);

    CHECK(body.dimension() == dim);

    // Origin should be inside (inside ball and x[0] - 0.5 <= 0)
    Point origin(dim);
    CHECK(body.is_in(origin) == -1);

    // Point inside ball but beyond hyperplane (x[0] = 0.75)
    // ball: 0.75^2 + 0^2 - 1 = -0.4375 <= 0, inside
    // hyperplane: 0.75 - 0.5 = 0.25 > 0, outside
    Point beyond_hyperplane(dim);
    beyond_hyperplane.set_coord(0, 0.75);
    CHECK(body.is_in(beyond_hyperplane) == 0);

    // Point outside ball
    Point far_outside(dim);
    far_outside.set_coord(0, 0.4);
    far_outside.set_coord(1, 1.5);
    CHECK(body.is_in(far_outside) == 0);
}

template <typename NT>
void test_line_positive_intersect() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef ConvexBody<Point> CB;

    unsigned int dim = 2;
    CB ball = generate_unit_ball<CB>(dim);

    // Intersect from origin along positive x axis (1, 0)
    // For unit ball: ||t * (1,0)||^2 - 1 = 0 => t = 1
    Point origin(dim);
    Point direction(dim);
    direction.set_coord(0, 1.0);

    auto result = ball.line_positive_intersect(origin, direction);
    NT t = result.first;

    CHECK(std::abs(t - 1.0) < 1e-4);
    // binary_search only searches in [0, 1], so for the unit ball
    // from the origin, the first constraint (ball) returns t ≈ 1.0
}

template <typename NT>
void test_binary_search() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef typename Point::FT FT;
    typedef ConvexBody<Point> CB;
    typedef std::function<NT(const Point&)> func;

    unsigned int dim = 2;

    // Linear constraint: x[0] - 0.3 <= 0
    func linear = [](const Point &x) {
        return x[0] - 0.3;
    };

    std::vector<func> gs{linear};
    std::vector< std::function<Point(const Point&)> > grads{
        [dim](const Point &x) {
            Point g(dim);
            g.set_coord(0, 1.0);
            return g;
        }
    };

    CB body(gs, grads, dim);

    // From origin along (1, 0): f(origin + t * (1,0)) = t - 0.3 = 0 => t = 0.3
    Point origin(dim);
    Point direction(dim);
    direction.set_coord(0, 1.0);

    NT t = body.binary_search(origin, direction, linear);
    CHECK(std::abs(t - 0.3) < 1e-4);
}

TEST_CASE("convex_body") {
    test_generate_unit_ball<double>();

    test_generate_unit_ball_intersect_hyperplane<double>();

    test_line_positive_intersect<double>();

    test_binary_search<double>();
}
