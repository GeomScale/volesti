// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis
// Copyright (c) 2024 Elias Tsigaridas

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef FEASIBLE_POINT_HPP
#define FEASIBLE_POINT_HPP

#include <tuple>

#include "preprocess/max_inscribed_ball.hpp"
#include "sampling/sphere.hpp" 
#include "convex_bodies/hpolytope.h"


// Using MT as to deal with both dense and sparse matrices
template <typename MT, typename VT>
VT compute_feasible_point(MT const& A, VT const& b)
{
    VT x;
    bool feasibility_only = true, converged;
    unsigned max_iters = 10000;
    // Compute a feasible point
    std::tie(x, std::ignore, converged) = max_inscribed_ball(A, b, max_iters, 1e-08, feasibility_only);
    if (!converged || ((A * x).array() > b.array()).any())
    {
        std::runtime_error("The computation of a feasible point failed.");
    }
    return x;
}

template <typename Point, typename Polytope, typename RNG>
std::pair<typename Polytope::VT,int>  compute_boundary_point(Polytope const& P, RNG& rng, typename Point::FT eps)
{
    using VT = typename Polytope::VT;
    using NT = typename Point::FT;
    const std::size_t m = P.num_of_hyperplanes();

    // Find the interior point
    VT r = compute_feasible_point(P.get_mat(), P.get_vec());

    // Random ray 
    const int dim = P.dimension();
    Point v_pt = GetDirection<Point>::apply(dim, rng);
    VT v = v_pt.getCoefficients();

    // First‐hit oracle
    VT Ar(m), Av(m);
    struct Params { NT inner_vi_ak; int facet_prev; };
    Params params;

    auto [lambda_min, facet] = P.line_first_positive_intersect(r, v, Ar, Av, params);

    // Checks
    if (!std::isfinite(lambda_min) || lambda_min <= eps || facet < 0)
        throw std::runtime_error("Failed to hit boundary!!!");

    // Compute boundry point + final check
    VT x = r + lambda_min * v;
    if ((P.get_mat() * x - P.get_vec()).maxCoeff() > eps)
        throw std::runtime_error("Boundary point violates constraints!!!");

    return std::pair<VT,int>(x, facet);    
}

#endif
