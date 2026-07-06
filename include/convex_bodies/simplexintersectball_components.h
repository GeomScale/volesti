#ifndef SIMPLEXINTERSECTBALL_COMPONENTS_H
#define SIMPLEXINTERSECTBALL_COMPONENTS_H

#include <cmath>
#include <queue>
#include <utility>
#include <vector>

#include <Eigen/Eigen>

/// Returns true if the segment [u, v] intersects the ball B(center, radius).
template <typename NT>
bool segment_intersects_ball(
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& u,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& v,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    VT direction = v - u;
    VT shifted = u - center;

    NT a = direction.dot(direction);
    NT b = NT(2) * shifted.dot(direction);
    NT c = shifted.dot(shifted) - radius * radius;

    if (a <= tol)
    {
        return c <= tol;
    }

    NT discriminant = b * b - NT(4) * a * c;

    if (discriminant < -tol)
    {
        return false;
    }

    if (discriminant < NT(0))
    {
        discriminant = NT(0);
    }

    NT sqrt_discriminant = std::sqrt(discriminant);

    NT t1 = (-b - sqrt_discriminant) / (NT(2) * a);
    NT t2 = (-b + sqrt_discriminant) / (NT(2) * a);

    return (t1 >= -tol && t1 <= NT(1) + tol) ||
           (t2 >= -tol && t2 <= NT(1) + tol);
}

/// Returns true if p lies strictly inside the ball B(center, radius).
template <typename NT>
bool point_is_inside_ball(
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& p,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    return (p - center).squaredNorm() < radius * radius - tol;
}

/// Returns a mask for simplex vertices that are outside the ball.
template <typename NT>
std::vector<int> active_vertices_outside_ball(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    int n = vertices.cols();
    std::vector<int> active(n, 0);

    for (int i = 0; i < n; ++i)
    {
        VT vertex = vertices.col(i);

        if (!point_is_inside_ball(vertex, center, radius, tol))
        {
            active[i] = 1;
        }
    }

    return active;
}

/// Builds the graph used to identify connected components of the
/// simplex-sphere intersection.
template <typename NT>
Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> build_simplex_ball_graph(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    int n = vertices.cols();

    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> adjacency =
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);

    std::vector<int> active =
        active_vertices_outside_ball(vertices, center, radius, tol);

    for (int i = 0; i < n; ++i)
    {
        if (!active[i])
        {
            continue;
        }

        VT vi = vertices.col(i);

        for (int j = i + 1; j < n; ++j)
        {
            if (!active[j])
            {
                continue;
            }

            VT vj = vertices.col(j);

            if (!segment_intersects_ball(vi, vj, center, radius, tol))
            {
                adjacency(i, j) = 1;
                adjacency(j, i) = 1;
            }
        }
    }

    return adjacency;
}

/// Finds connected components of adjacency restricted to active vertices.
inline std::vector<std::vector<int>> connected_components_from_graph(
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> const& adjacency,
    std::vector<int> const& active)
{
    int n = adjacency.rows();

    std::vector<int> visited(n, 0);
    std::vector<std::vector<int>> components;

    for (int start = 0; start < n; ++start)
    {
        if (visited[start] || !active[start])
        {
            continue;
        }

        std::vector<int> component;
        std::queue<int> queue;

        visited[start] = 1;
        queue.push(start);

        while (!queue.empty())
        {
            int current = queue.front();
            queue.pop();

            component.push_back(current);

            for (int next = 0; next < n; ++next)
            {
                if (!visited[next] && active[next] && adjacency(current, next) != 0)
                {
                    visited[next] = 1;
                    queue.push(next);
                }
            }
        }

        components.push_back(component);
    }

    return components;
}

/// Finds connected components of adjacency, treating all vertices as active.
inline std::vector<std::vector<int>> connected_components_from_graph(
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> const& adjacency)
{
    int n = adjacency.rows();

    std::vector<int> active(n, 1);

    return connected_components_from_graph(adjacency, active);
}

/// Finds connected components of the simplex-sphere intersection.
template <typename NT>
std::vector<std::vector<int>> find_simplex_ball_components(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> adjacency =
        build_simplex_ball_graph(vertices, center, radius, tol);

    std::vector<int> active =
        active_vertices_outside_ball(vertices, center, radius, tol);

    return connected_components_from_graph(adjacency, active);
}

/// Intersects the ray from an interior point to a vertex with the sphere boundary.
template <typename NT>
std::pair<bool, Eigen::Matrix<NT, Eigen::Dynamic, 1>>
ray_sphere_intersection_from_interior(
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& interior_point,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& vertex,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    VT direction = vertex - interior_point;
    VT shifted = interior_point - center;

    NT a = direction.dot(direction);
    NT b = NT(2) * shifted.dot(direction);
    NT c = shifted.dot(shifted) - radius * radius;

    if (a <= tol)
    {
        VT empty(center.rows());
        empty.setZero();

        return std::make_pair(false, empty);
    }

    NT discriminant = b * b - NT(4) * a * c;

    if (discriminant < -tol)
    {
        VT empty(center.rows());
        empty.setZero();

        return std::make_pair(false, empty);
    }

    if (discriminant < NT(0))
    {
        discriminant = NT(0);
    }

    NT sqrt_discriminant = std::sqrt(discriminant);

    NT t1 = (-b - sqrt_discriminant) / (NT(2) * a);
    NT t2 = (-b + sqrt_discriminant) / (NT(2) * a);

    NT t = t1;

    if (t < -tol || t > NT(1) + tol)
    {
        t = t2;
    }

    if (t < -tol || t > NT(1) + tol)
    {
        VT empty(center.rows());
        empty.setZero();

        return std::make_pair(false, empty);
    }

    VT candidate = interior_point + t * direction;

    return std::make_pair(true, candidate);
}

/// Returns true if p satisfies A * p <= b.
template <typename NT>
bool point_satisfies_halfspaces(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& p,
    NT tol = NT(1e-10))
{
    for (int i = 0; i < A.rows(); ++i)
    {
        if (A.row(i).dot(p) > b(i) + tol)
        {
            return false;
        }
    }

    return true;
}

/// Finds a starting point for one component by intersecting rays from an
/// interior point to the component vertices with the sphere boundary.
template <typename NT>
std::pair<bool, Eigen::Matrix<NT, Eigen::Dynamic, 1>>
find_starting_point_for_component(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    std::vector<int> const& component,
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& interior_point,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    for (int vertex_index : component)
    {
        VT vertex = vertices.col(vertex_index);

        std::pair<bool, VT> intersection =
            ray_sphere_intersection_from_interior(
                interior_point, vertex, center, radius, tol);

        if (!intersection.first)
        {
            continue;
        }

        VT candidate = intersection.second;

        bool on_sphere =
            std::abs((candidate - center).norm() - radius) <= NT(100) * tol;

        if (on_sphere && point_satisfies_halfspaces(A, b, candidate, tol))
        {
            return std::make_pair(true, candidate);
        }
    }

    VT empty(center.rows());
    empty.setZero();

    return std::make_pair(false, empty);
}

/// Finds one starting point for each connected component.
template <typename NT>
std::vector<Eigen::Matrix<NT, Eigen::Dynamic, 1>>
find_starting_points_for_components(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    std::vector<std::vector<int>> const& components,
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& interior_point,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    std::vector<VT> starting_points;

    for (std::vector<int> const& component : components)
    {
        std::pair<bool, VT> result =
            find_starting_point_for_component(
                vertices, component, A, b, interior_point, center, radius, tol);

        if (result.first)
        {
            starting_points.push_back(result.second);
        }
    }

    return starting_points;
}

/// Finds connected components and corresponding starting points.
template <typename NT>
std::pair<
    std::vector<std::vector<int>>,
    std::vector<Eigen::Matrix<NT, Eigen::Dynamic, 1>>>
find_simplex_ball_components_and_starting_points(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& interior_point,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    std::vector<std::vector<int>> components =
        find_simplex_ball_components(vertices, center, radius, tol);

    std::vector<Eigen::Matrix<NT, Eigen::Dynamic, 1>> starting_points =
        find_starting_points_for_components(
            vertices, components, A, b, interior_point, center, radius, tol);

    return std::make_pair(components, starting_points);
}

#endif