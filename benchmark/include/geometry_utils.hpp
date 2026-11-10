#pragma once

#include "core_types.hpp"
#include <cmath>

// Useful function to rotate polytope
template <typename HPOLYTOPE>
HPOLYTOPE rotate_all_dims(const HPOLYTOPE& P, typename HPOLYTOPE::NT angle)
{
    using NT = typename HPOLYTOPE::NT;
    int dim = P.dimension();

    // Build global rotation matrix
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> R =
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>::Identity(dim, dim);

    NT c = std::cos(angle);
    NT s = std::sin(angle);

    // Apply rotation in each adjacent coordinate plane
    for (int k = 0; k < dim - 1; ++k) {
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> Rk =
            Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>::Identity(dim, dim);

        Rk(k,   k)   =  c;
        Rk(k,   k+1) = -s;
        Rk(k+1, k)   =  s;
        Rk(k+1, k+1) =  c;

        R = R * Rk;   // compose rotations
    }

    // Extract A and b
    auto A = P.get_mat();
    auto b = P.get_vec();

    // Apply A' = A R
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> A_rot = A * R;

    return HPOLYTOPE(dim, A_rot, b);
}