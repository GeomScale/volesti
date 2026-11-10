// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONSTRAINED_POLYTOPE_PRESOLVE_HPP
#define CONSTRAINED_POLYTOPE_PRESOLVE_HPP

#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <limits>
#include <Eigen/Sparse>
#include "cholmod.h"
#include "SuiteSparseQR.hpp"
#include "convex_bodies/constrained_polytope.hpp"
#include "preprocess/constrained_polytope/dependent_rows.hpp"


// Result of finding the variables pinned by equalities.
template <typename Point>
struct PinnedPresolveResult {
    // The polytope with the bounds freed.
    ConstrainedPolytope<Point> polytope;

    // The pinned variables and their values.
    std::vector<unsigned> pinned;

    // The rank of A_eq.
    unsigned rank = 0;

    // False if the presolve step failed.
    bool valid = false;
};

// Frees the bounds of every variable that A_eq x = b_eq pins alone.
// @tparam Point the point type
// @param P the polytope
// @param tol the rank zero tolerance
// @return the polytope with the bounds of pinned variables freed
template <typename Point>
PinnedPresolveResult<Point> presolve_pinned_vars(ConstrainedPolytope<Point> const& P,
                                                 double tol = 1e-7)
{
    using NT = typename ConstrainedPolytope<Point>::NT;
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;

    PinnedPresolveResult<Point> res;

    // Quick check
    if (P.getNumEqualities() == 0) {
        res.polytope = P;
        res.valid = true;
        return res;
    }

    unsigned const d = P.getDimension();
    MT const& A_eq = P.getEqualities();
    VT const& b_eq = P.getEqualityRHS();
    Eigen::SparseMatrix<double> A = A_eq.template cast<double>();

    cholmod_common cc;
    cholmod_l_start(&cc);

    Eigen::MatrixXd N;
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(d);

    try {
        N = dependent_rows_removal_util::nullspace_basis(A, tol, &cc, res.rank);

        std::vector<unsigned> keep = dependent_rows_removal_util::independent_rows(A, tol, &cc);

        if (!keep.empty()) {
            std::vector<Eigen::Triplet<double>> triplets;
            Eigen::VectorXd b_keep((Eigen::Index)keep.size());

            for (std::size_t r = 0; r < keep.size(); ++r) {
                for (typename MT::InnerIterator it(A_eq, keep[r]); it; ++it)
                    triplets.emplace_back((int)r, (int)it.col(), (double)it.value());
                b_keep((Eigen::Index)r) = (double)b_eq(keep[r]);
            }

            Eigen::SparseMatrix<double> A_keep((Eigen::Index)keep.size(), d);
            A_keep.setFromTriplets(triplets.begin(), triplets.end());

            x0 = dependent_rows_removal_util::min_norm_solution(A_keep, b_keep, tol, &cc);
        }
        cholmod_l_finish(&cc);
    } catch (...) {
        cholmod_l_finish(&cc);
        throw;
    }

    NT const INF = std::numeric_limits<NT>::infinity();
    VT lb = P.getLowerBounds();
    VT ub = P.getUpperBounds();

    for (unsigned j = 0; j < d; ++j) {
        if (N.cols() > 0 && N.row(j).cwiseAbs().maxCoeff() >= tol) continue;

        if (x0(j) < (double)lb(j)-tol || x0(j) > (double)ub(j)+tol) return res;

        res.pinned.push_back(j);
        lb(j) = -INF;
        ub(j) = INF;
    }

    res.polytope = ConstrainedPolytope<Point>(d, A_eq, b_eq, P.getInequalities(), P.getInequalityRHS(),
                                              lb, ub);
    
    res.valid = true;
    return res;
}
#endif