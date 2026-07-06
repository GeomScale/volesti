#include "doctest.h"

#include <Eigen/Eigen>

#include "convex_bodies/simplexintersectball_components.h"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

TEST_CASE("simplexball_components_segment_intersects_ball")
{
    VT center(3);
    center << 0, 0, 0;

    VT u(3);
    VT v(3);

    // Segment crosses the unit ball.
    u << -2, 0, 0;
    v << 2, 0, 0;
    CHECK(segment_intersects_ball(u, v, center, NT(1)));

    // Segment stays outside the unit ball.
    u << 2, 2, 0;
    v << 3, 2, 0;
    CHECK_FALSE(segment_intersects_ball(u, v, center, NT(1)));

    // Segment is tangent to the unit ball.
    u << -1, 1, 0;
    v << 1, 1, 0;
    CHECK(segment_intersects_ball(u, v, center, NT(1)));
}

TEST_CASE("simplexball_components_tetrahedron_finds_two_components")
{
    VT center(3);
    center << 0, 0, 0;

    // Non-degenerate tetrahedron in R^3.
    // Vertices 0 and 1 form one component.
    // Vertices 2 and 3 form the other component.
    // Every edge between the two pairs intersects the unit ball.
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vertices(3, 4);
    vertices << -2, -2,  2,  2,
                -0.2, 0.2, -0.2, 0.2,
                -0.2, 0.2,  0.2, -0.2;

    std::vector<std::vector<int>> components =
        find_simplex_ball_components(vertices, center, NT(1));

    REQUIRE(components.size() == 2);

    CHECK(components[0].size() == 2);
    CHECK(components[0][0] == 0);
    CHECK(components[0][1] == 1);

    CHECK(components[1].size() == 2);
    CHECK(components[1][0] == 2);
    CHECK(components[1][1] == 3);
}

TEST_CASE("simplexball_components_tetrahedron_starting_points")
{
    VT center(3);
    center << 0, 0, 0;

    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vertices(3, 4);
    vertices << -2, -2,  2,  2,
                -0.2, 0.2, -0.2, 0.2,
                -0.2, 0.2,  0.2, -0.2;

    // H-representation of the tetrahedron above.
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> A(4, 3);
    A <<  1,  10,  10,
          1, -10, -10,
         -1,  10, -10,
         -1, -10,  10;

    VT b(4);
    b << 2, 2, 2, 2;

    VT interior_point(3);
    interior_point << 0, 0, 0;

    std::pair<std::vector<std::vector<int>>, std::vector<VT>> result =
        find_simplex_ball_components_and_starting_points(
            vertices, A, b, interior_point, center, NT(1));

    REQUIRE(result.first.size() == 2);
    REQUIRE(result.second.size() == 2);

    for (VT const& p : result.second)
    {
        CHECK(p.rows() == 3);
        CHECK(p.norm() == doctest::Approx(1.0));
        CHECK(point_satisfies_halfspaces(A, b, p));
    }
}

TEST_CASE("simplexball_components_tetrahedron_filters_interior_vertex")
{
    VT center(3);
    center << 0, 0, 0;

    // Non-degenerate tetrahedron with one vertex inside the unit ball.
    // Vertex 3 is inside and should be ignored.
    // Vertex 0 remains an isolated active component.
    // Vertices 1 and 2 remain connected.
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vertices(3, 4);
    vertices << -2,  2,  2,  0,
                 0,  0,  0.5, 0.1,
                 0,  0,  0.5, 0;

    std::vector<std::vector<int>> components =
        find_simplex_ball_components(vertices, center, NT(1));

    REQUIRE(components.size() == 2);

    CHECK(components[0].size() == 1);
    CHECK(components[0][0] == 0);

    CHECK(components[1].size() == 2);
    CHECK(components[1][0] == 1);
    CHECK(components[1][1] == 2);
}