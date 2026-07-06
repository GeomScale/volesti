#ifndef SIMPLEXINTERSECTBALL_COMPONENTS_H
#define SIMPLEXINTERSECTBALL_COMPONENTS_H

#include <cmath>
#include <vector>
#include <queue>
#include <utility>
#include <Eigen/Eigen>

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

template <typename NT>
bool point_is_inside_ball(
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& p,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    return (p - center).squaredNorm() < radius * radius - tol;
}

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
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Ones(n, n);

    for (int i = 0; i < n; ++i)
    {
        adjacency(i, i) = 0;
    }

    for (int i = 0; i < n; ++i)
    {
        VT vi = vertices.col(i);

        if (point_is_inside_ball(vi, center, radius, tol))
        {
            adjacency.row(i).setZero();
            adjacency.col(i).setZero();
            continue;
        }

        for (int j = i + 1; j < n; ++j)
        {
            VT vj = vertices.col(j);

            if (point_is_inside_ball(vj, center, radius, tol) ||
                segment_intersects_ball(vi, vj, center, radius, tol))
            {
                adjacency(i, j) = 0;
                adjacency(j, i) = 0;
            }
        }
    }

    return adjacency;
}

inline std::vector<std::vector<int>> connected_components_from_graph(
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> const& adjacency)
{
    int n = adjacency.rows();

    std::vector<int> visited(n, 0);
    std::vector<std::vector<int>> components;

    for (int start = 0; start < n; ++start)
    {
        if (visited[start])
        {
            continue;
        }

        bool isolated_removed_vertex = true;
        for (int j = 0; j < n; ++j)
        {
            if (adjacency(start, j) != 0 || adjacency(j, start) != 0)
            {
                isolated_removed_vertex = false;
                break;
            }
        }

        if (isolated_removed_vertex)
        {
            visited[start] = 1;
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
                if (!visited[next] && adjacency(current, next) != 0)
                {
                    visited[next] = 1;
                    queue.push(next);
                }
            }
        }

        if (!component.empty())
        {
            components.push_back(component);
        }
    }

    return components;
}

template <typename NT>
std::vector<std::vector<int>> find_simplex_ball_components(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> adjacency =
        build_simplex_ball_graph(vertices, center, radius, tol);

    return connected_components_from_graph(adjacency);
}

template <typename NT>
Eigen::Matrix<NT, Eigen::Dynamic, 1> radial_starting_point_from_vertex(
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& vertex,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    VT direction = vertex - center;
    NT norm = direction.norm();

    if (norm <= tol)
    {
        return center;
    }

    return center + radius * direction / norm;
}

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

template <typename NT>
std::pair<bool, Eigen::Matrix<NT, Eigen::Dynamic, 1>>
find_starting_point_for_component(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    std::vector<int> const& component,
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    for (int vertex_index : component)
    {
        VT vertex = vertices.col(vertex_index);
        VT candidate = radial_starting_point_from_vertex(vertex, center, radius, tol);

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

template <typename NT>
std::vector<Eigen::Matrix<NT, Eigen::Dynamic, 1>>
find_starting_points_for_components(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    std::vector<std::vector<int>> const& components,
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    std::vector<VT> starting_points;

    for (std::vector<int> const& component : components)
    {
        std::pair<bool, VT> result =
            find_starting_point_for_component(vertices, component, A, b, center, radius, tol);

        if (result.first)
        {
            starting_points.push_back(result.second);
        }
    }

    return starting_points;
}

template <typename NT>
std::pair<
    std::vector<std::vector<int>>,
    std::vector<Eigen::Matrix<NT, Eigen::Dynamic, 1>>>
find_simplex_ball_components_and_starting_points(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& vertices,
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    std::vector<std::vector<int>> components =
        find_simplex_ball_components(vertices, center, radius, tol);

    std::vector<Eigen::Matrix<NT, Eigen::Dynamic, 1>> starting_points =
        find_starting_points_for_components(vertices, components, A, b, center, radius, tol);

    return std::make_pair(components, starting_points);
}

#endif