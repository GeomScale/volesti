// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RAY_SHOOT_HPP
#define RAY_SHOOT_HPP

#include <cmath>
#include <vector>
#include <limits>
#include "convex_bodies/constrained_polytope.hpp"

// The result of a ray shoot.
struct RayShoot {
    // The index of the first inequality hit.
    unsigned row = 0;

    // The t for which the ray hits the facet.
    double t = std::numeric_limits<double>::infinity();

    // True if a row was hit.
    bool hit = false;
};

// Shoots the ray z + t*r, t>=0 and returns the first facet (row of A_in)
// it hits.
//
// Only rows in the candidate set J are tested, this is an optimization of Clarkson.
//
//      a_i^T r > tol
//      t = (b_in(i) - a_i^T z) / (a_i^T r) >= 0
//
// @tparam Point the point type
// @param P the polytope
// @param J the indices of the candidate rows to test
// @param z the ray origin, must be strictly interior
// @param r the ray direction
// @param tol directions smaller than this are treated as zero
// @return the ray shoot result
template <typename Point>
RayShoot ray_shoot(ConstrainedPolytope<Point> const& P,
                   std::vector<unsigned> const& J,
                   typename ConstrainedPolytope<Point>::VT const& z,
                   typename ConstrainedPolytope<Point>::VT const& r,
                   double tol = 1e-9)
{
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;
    
    MT const& A_in = P.getInequalities();
    VT const& b_in = P.getInequalityRHS();

    RayShoot rs;

    for (unsigned i : J) {
        // How fast the ray moves towards i
        double ar = (double)A_in.row(i).dot(r);
        if (ar <= tol) continue;

        // Where the ray starts from
        double az = (double)A_in.row(i).dot(z);
        double t = ((double)b_in(i)-az)/ar;

        if (t < 0.0) continue;

        if (!rs.hit || rs.t > t) {
            rs.t = t;
            rs.row = i;
            rs.hit = true;
        }

    }

    return rs;
}
#endif