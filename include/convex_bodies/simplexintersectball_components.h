#ifndef SIMPLEXINTERSECTBALL_COMPONENTS_H
#define SIMPLEXINTERSECTBALL_COMPONENTS_H

#include <cmath>
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

#endif