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