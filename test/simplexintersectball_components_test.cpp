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

TEST_CASE("simplexball_components_connected_components_from_graph")
{
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> graph(5, 5);
    graph.setZero();

    // Component 1: 0 -- 1
    graph(0, 1) = 1;
    graph(1, 0) = 1;

    // Component 2: 2 -- 3
    graph(2, 3) = 1;
    graph(3, 2) = 1;

    // Vertex 4 is isolated and should be ignored.

    std::vector<std::vector<int>> components =
        connected_components_from_graph(graph);

    CHECK(components.size() == 2);

    CHECK(components[0].size() == 2);
    CHECK(components[1].size() == 2);

    CHECK(components[0][0] == 0);
    CHECK(components[0][1] == 1);

    CHECK(components[1][0] == 2);
    CHECK(components[1][1] == 3);
}

TEST_CASE("simplexball_components_find_components_from_vertices")
{
    VT center(2);
    center << 0, 0;

    // Same 4 vertices as before:
    // left, right, top, bottom around the unit ball.
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vertices(2, 4);
    vertices << -2,  2,  0,  0,
                 0,  0,  2, -2;

    std::vector<std::vector<int>> components =
        find_simplex_ball_components(vertices, center, NT(1));

    CHECK(components.size() == 1);

    CHECK(components[0].size() == 4);
    CHECK(components[0][0] == 0);
    CHECK(components[0][1] == 2);
    CHECK(components[0][2] == 3);
    CHECK(components[0][3] == 1);
}

TEST_CASE("simplexball_components_radial_starting_point_from_vertex")
{
    VT center(2);
    center << 0, 0;

    VT vertex(2);
    vertex << 2, 0;

    VT p = radial_starting_point_from_vertex(vertex, center, NT(1));

    CHECK(p.rows() == 2);
    CHECK(p(0) == doctest::Approx(1.0));
    CHECK(p(1) == doctest::Approx(0.0));
    CHECK(p.norm() == doctest::Approx(1.0));

    vertex << 0, -3;

    p = radial_starting_point_from_vertex(vertex, center, NT(1));

    CHECK(p(0) == doctest::Approx(0.0));
    CHECK(p(1) == doctest::Approx(-1.0));
    CHECK(p.norm() == doctest::Approx(1.0));
}

TEST_CASE("simplexball_components_find_starting_point_for_component")
{
    VT center(2);
    center << 0, 0;

    // Triangle:
    // x >= 0, y >= 0, x + y <= 2
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> A(3, 2);
    A << -1,  0,
          0, -1,
          1,  1;

    VT b(3);
    b << 0, 0, 2;

    // Vertices stored column-wise: (0,0), (2,0), (0,2)
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vertices(2, 3);
    vertices << 0, 2, 0,
                0, 0, 2;

    std::vector<int> component;
    component.push_back(1); // vertex (2,0)
    component.push_back(2); // vertex (0,2)

    std::pair<bool, VT> result =
        find_starting_point_for_component(vertices, component, A, b, center, NT(1));

    CHECK(result.first);

    VT p = result.second;

    CHECK(p.norm() == doctest::Approx(1.0));
    CHECK(point_satisfies_halfspaces(A, b, p));
}

TEST_CASE("simplexball_components_find_starting_points_for_components")
{
    VT center(2);
    center << 0, 0;

    // Triangle:
    // x >= 0, y >= 0, x + y <= 2
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> A(3, 2);
    A << -1,  0,
          0, -1,
          1,  1;

    VT b(3);
    b << 0, 0, 2;

    // Vertices stored column-wise: (0,0), (2,0), (0,2)
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vertices(2, 3);
    vertices << 0, 2, 0,
                0, 0, 2;

    std::vector<std::vector<int>> components;
    components.push_back({1, 2});

    std::vector<VT> starting_points =
        find_starting_points_for_components(vertices, components, A, b, center, NT(1));

    CHECK(starting_points.size() == 1);

    VT p = starting_points[0];

    CHECK(p.norm() == doctest::Approx(1.0));
    CHECK(point_satisfies_halfspaces(A, b, p));
}