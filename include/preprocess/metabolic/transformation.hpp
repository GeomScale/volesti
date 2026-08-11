// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_TRANSFORMATION_HPP
#define METABOLIC_TRANSFORMATION_HPP

#include "convex_bodies/metabolic_polytope.hpp"
#include "convex_bodies/hpolytope.h"
#include "preprocess/full_dimensional_polytope.hpp"

// Transforms a MetabolicPolytope into a full-dimensional HPolytope.
// @tparam Point the point type of the polytope
// @param P the input MetabolicPolytope
// @return a tuple (HPolytope, shift, N):
//         - HPolytope : the full-dimensional polytope
//         - shift     : a particular solution of the A_eq x = b_eq
//         - N         : a basis of the nullspace of A_eq
template <typename Point>
std::tuple<HPolytope<Point>,
            typename MetabolicPolytope<Point>::VT,
            typename Eigen::Matrix<typename MetabolicPolytope<Point>::NT, 
                                    Eigen::Dynamic, Eigen::Dynamic>>
transform(MetabolicPolytope<Point> const& P) {
    typedef typename MetabolicPolytope<Point>::MT MT;
    typedef typename MetabolicPolytope<Point>::VT VT;
    typedef typename MetabolicPolytope<Point>::NT NT;
    typedef typename Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DenseMT;
    typedef typename Eigen::SparseMatrix<NT, Eigen::ColMajor> MTColMajor;

    unsigned d = P.getDimension();

    MTColMajor A_eq = P.getEqualities();
    const VT& b_eq = P.getEqualityBounds();
    const VT& b_l = P.getLowerBounds();
    const VT& b_u = P.getUpperBounds();

    unsigned rows = P.getNumFiniteBounds();

    DenseMT A(rows, d);
    A.setZero();
    VT b(rows);

    unsigned l = 0;
    for (unsigned k = 0; k < d; ++k) {
        if (!std::isinf(b_l(k))) {
            A(l, k) = NT(-1.0);
            b(l) = -b_l(k);
            ++l;
        }
        if (!std::isinf(b_u(k))) {
            A(l, k) = NT(1.0);
            b(l) = b_u(k);
            ++l;
        }
    }

    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope<NT, MTColMajor, DenseMT, VT>(
        A_eq, 
        b_eq, 
        A, 
        b
    );

    HPolytope<Point> HP(A_full.cols(), A_full, b_full);
    return std::make_tuple(HP, shift, N);
}
#endif