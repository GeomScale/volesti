#ifndef SIMPLEXINTERSECTBALL_COMPONENTS_H
#define SIMPLEXINTERSECTBALL_COMPONENTS_H

#include <cmath>
#include <cstdlib>
#include <queue>
#include <utility>
#include <vector>

#include <Eigen/Eigen>

#undef Realloc
#undef Free
#include "lp_lib.h"


/// Solves ||point + t * direction - center||^2 = radius^2.
/// Returns false if the line does not intersect the sphere.
/// Precondition: direction.dot(direction) > tol.
template <typename NT>
bool solve_ball_line_roots(
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& point,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& direction,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius,
    NT tol,
    NT& tmin,
    NT& tmax)
{
    Eigen::Matrix<NT, Eigen::Dynamic, 1> shifted = point - center;

    NT a = direction.dot(direction);
    NT b = NT(2) * shifted.dot(direction);
    NT c = shifted.dot(shifted) - radius * radius;

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

    tmin = (-b - sqrt_discriminant) / (NT(2) * a);
    tmax = (-b + sqrt_discriminant) / (NT(2) * a);

    return true;
}

/// Tests whether a segment intersects a Euclidean ball.
template <typename NT>
bool segment_intersects_ball(
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& u,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& v,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    NT tol = NT(1e-10))
{
    Eigen::Matrix<NT, Eigen::Dynamic, 1> direction = v - u;

    if (direction.dot(direction) <= tol)
    {
        return (u - center).squaredNorm() <= radius * radius + tol;
    }

    NT tmin;
    NT tmax;

    if (!solve_ball_line_roots(u, direction, center, radius, tol, tmin, tmax))
    {
        return false;
    }

    return (tmin >= -tol && tmin <= NT(1) + tol) ||
           (tmax >= -tol && tmax <= NT(1) + tol);
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

    VT empty = VT::Zero(center.rows());
    VT direction = vertex - interior_point;

    if (direction.dot(direction) <= tol)
    {
        return std::make_pair(false, empty);
    }

    NT tmin;
    NT tmax;

    if (!solve_ball_line_roots(
            interior_point, direction, center, radius, tol, tmin, tmax))
    {
        return std::make_pair(false, empty);
    }

    NT t = tmin;

    if (t < -tol || t > NT(1) + tol)
    {
        t = tmax;
    }

    if (t < -tol || t > NT(1) + tol)
    {
        return std::make_pair(false, empty);
    }

    VT candidate = interior_point + t * direction;

    return std::make_pair(true, candidate);
}

/// Computes an approximate Chebyshev center of {x : A x <= b} intersected
/// with B(center, radius), using cutting-plane linearization of the ball.
template <typename NT>
std::pair<bool, Eigen::Matrix<NT, Eigen::Dynamic, 1>>
chebyshev_center_intersect_ball(
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> const& A,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& b,
    Eigen::Matrix<NT, Eigen::Dynamic, 1> const& center,
    NT radius = NT(1),
    unsigned int max_iterations = 50,
    NT tol = NT(1e-8))
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    int m = A.rows();
    int d = A.cols();
    int ncols = d + 1;

    std::vector<VT> ball_cuts;

    for (int i = 0; i < d; ++i)
    {
        VT positive = VT::Zero(d);
        positive(i) = NT(1);
        ball_cuts.push_back(positive);

        VT negative = VT::Zero(d);
        negative(i) = NT(-1);
        ball_cuts.push_back(negative);
    }

    VT solution = center;

    for (unsigned int iteration = 0; iteration < max_iterations; ++iteration)
    {
        lprec* lp = make_lp(0, ncols);

        if (lp == NULL)
        {
            return std::make_pair(false, solution);
        }

        REAL infinite = get_infinite(lp);

        for (int j = 0; j < d; ++j)
        {
            set_bounds(lp, j + 1, -infinite, infinite);
        }

        set_bounds(lp, d + 1, 0.0, infinite);
        set_add_rowmode(lp, TRUE);

        std::vector<int> colno(ncols);
        std::vector<REAL> row(ncols);

        for (int j = 0; j < ncols; ++j)
        {
            colno[j] = j + 1;
        }

        for (int i = 0; i < m; ++i)
        {
            NT normal_norm = A.row(i).norm();

            for (int j = 0; j < d; ++j)
            {
                row[j] = A(i, j);
            }

            row[d] = normal_norm;

            if (!add_constraintex(lp, ncols, row.data(), colno.data(), LE, b(i)))
            {
                delete_lp(lp);
                return std::make_pair(false, solution);
            }
        }

        for (VT const& cut : ball_cuts)
        {
            for (int j = 0; j < d; ++j)
            {
                row[j] = cut(j);
            }

            row[d] = NT(1);

            NT rhs = radius + cut.dot(center);

            if (!add_constraintex(lp, ncols, row.data(), colno.data(), LE, rhs))
            {
                delete_lp(lp);
                return std::make_pair(false, solution);
            }
        }

        set_add_rowmode(lp, FALSE);

        for (int j = 0; j < d; ++j)
        {
            row[j] = NT(0);
        }

        row[d] = NT(1);

        if (!set_obj_fnex(lp, ncols, row.data(), colno.data()))
        {
            delete_lp(lp);
            return std::make_pair(false, solution);
        }

        set_maxim(lp);
        set_verbose(lp, NEUTRAL);

        if (solve(lp) != OPTIMAL)
        {
            delete_lp(lp);
            return std::make_pair(false, solution);
        }

        get_variables(lp, row.data());

        for (int j = 0; j < d; ++j)
        {
            solution(j) = row[j];
        }

        NT inner_radius = row[d];

        delete_lp(lp);

        VT shifted = solution - center;
        NT distance = shifted.norm();

        if (distance + inner_radius <= radius + tol)
        {
            return std::make_pair(true, solution);
        }

        if (distance <= tol)
        {
            return std::make_pair(false, solution);
        }

        ball_cuts.push_back(shifted / distance);
    }

    return std::make_pair(false, solution);
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

        if (point_satisfies_halfspaces(A, b, candidate, tol))
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

/// Finds connected components and starting points, computing an approximate
/// Chebyshev center of the simplex-ball intersection as the interior point.
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
    std::pair<bool, Eigen::Matrix<NT, Eigen::Dynamic, 1>> chebyshev_result =
        chebyshev_center_intersect_ball(A, b, center, radius, 50, tol);

    Eigen::Matrix<NT, Eigen::Dynamic, 1> interior_point =
        chebyshev_result.first ? chebyshev_result.second : center;

    return find_simplex_ball_components_and_starting_points(
        vertices, A, b, interior_point, center, radius, tol);
}

#endif