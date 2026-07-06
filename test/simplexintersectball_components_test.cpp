#include "doctest.h"

#include <Eigen/Eigen>

#include "convex_bodies/simplexintersectball_components.h"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

TEST_CASE("simplexball_components_segment_intersects_ball")
{
    VT center(2);
    center << 0, 0;

    VT u(2);
    VT v(2);

    // Segment crosses the unit ball.
    u << -2, 0;
    v << 2, 0;
    CHECK(segment_intersects_ball(u, v, center, NT(1)));

    // Segment stays outside the unit ball.
    u << 2, 2;
    v << 3, 2;
    CHECK_FALSE(segment_intersects_ball(u, v, center, NT(1)));

    // Segment is tangent to the unit ball.
    u << -1, 1;
    v << 1, 1;
    CHECK(segment_intersects_ball(u, v, center, NT(1)));

    // Segment starts inside the unit ball.
    u << 0, 0;
    v << 2, 0;
    CHECK(segment_intersects_ball(u, v, center, NT(1)));
}

TEST_CASE("simplexball_components_point_is_inside_ball")
{
    VT center(2);
    center << 0, 0;

    VT p(2);

    p << 0.5, 0.0;
    CHECK(point_is_inside_ball(p, center, NT(1)));

    p << 1.0, 0.0;
    CHECK_FALSE(point_is_inside_ball(p, center, NT(1)));

    p << 1.5, 0.0;
    CHECK_FALSE(point_is_inside_ball(p, center, NT(1)));
}

TEST_CASE("simplexball_components_build_graph")
{
    VT center(2);
    center << 0, 0;

    // 4 vertices around the unit ball.
    // Edges through the ball should be removed.
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vertices(2, 4);
    vertices << -2,  2,  0,  0,
                 0,  0,  2, -2;

    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> graph =
        build_simplex_ball_graph(vertices, center, NT(1));

    CHECK(graph.rows() == 4);
    CHECK(graph.cols() == 4);

    // No self-loops.
    CHECK(graph(0, 0) == 0);
    CHECK(graph(1, 1) == 0);
    CHECK(graph(2, 2) == 0);
    CHECK(graph(3, 3) == 0);

    // Segment from (-2,0) to (2,0) crosses the ball.
    CHECK(graph(0, 1) == 0);
    CHECK(graph(1, 0) == 0);

    // Segment from (0,2) to (0,-2) crosses the ball.
    CHECK(graph(2, 3) == 0);
    CHECK(graph(3, 2) == 0);

    // Segment from (-2,0) to (0,2) is tangent/outside boundary-connected.
    CHECK(graph(0, 2) == 1);
    CHECK(graph(2, 0) == 1);
}