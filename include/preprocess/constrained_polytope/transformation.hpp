// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONSTRAINED_POLYTOPE_TRANSFORMATION_HPP
#define CONSTRAINED_POLYTOPE_TRANSFORMATION_HPP

#include <cmath>
#include <tuple>
#include "convex_bodies/constrained_polytope.hpp"
#include "convex_bodies/hpolytope.h"
#include "preprocess/full_dimensional_polytope.hpp"

// Transforms a ConstrainedPolytope into a full-dimensional HPolytope.
// @tparam Point the point type of the polytope
// @param P the input ConstrainedPolytope
// @return a tuple (HPolytope, shift, N):
//         - HPolytope : the full-dimensional polytope
//         - shift     : a particular solution of the A_eq x = b_eq
//         - N         : a basis of the nullspace of A_eq
template <typename Point>
std::tuple<HPolytope<Point>,
            typename ConstrainedPolytope<Point>::VT,
            typename Eigen::Matrix<typename ConstrainedPolytope<Point>::NT, 
                                    Eigen::Dynamic, Eigen::Dynamic>>
transform(ConstrainedPolytope<Point> const& P) {
    typedef typename ConstrainedPolytope<Point>::MT MT;
    typedef typename ConstrainedPolytope<Point>::VT VT;
    typedef typename ConstrainedPolytope<Point>::NT NT;
    typedef typename Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DenseMT;
    typedef typename Eigen::SparseMatrix<NT, Eigen::ColMajor> MTColMajor;

    unsigned d = P.getDimension();

    MTColMajor A_eq = P.getEqualities();
    const VT& b_eq = P.getEqualityRHS();
    const MT& A_in = P.getInequalities();
    const VT& b_in = P.getInequalityRHS();
    const VT& b_l = P.getLowerBounds();
    const VT& b_u = P.getUpperBounds();
    unsigned const m_in = P.getNumInequalities();
    unsigned rows = m_in+P.getNumFiniteBounds();

    DenseMT A = DenseMT::Zero(rows, d);
    VT b(rows);

    // Adds the rows of A_in
    for (unsigned i = 0; i < m_in; ++i) {
        for (typename MT::InnerIterator it(A_in, i); it; ++it)
            A(i, it.col()) = it.value();
        b(i) = b_in(i);
    }

    // Adds the finite bounds
    unsigned l = m_in;
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

    // Without equalities the polytope is full dimensional
    if (A_eq.rows() == 0) {
        HPolytope<Point> HP(d, A, b);
        return std::make_tuple(HP, VT(VT::Zero(d)), DenseMT(DenseMT::Identity(d,d)));
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