// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HIGHS_POLYTOPE_HPP
#define HIGHS_POLYTOPE_HPP

#include <cmath>
#include <vector>
#include <limits>
#include "Highs.h"
#include "convex_bodies/constrained_polytope.hpp"

// Configures a HiGHS instance for use in simplification.
// @param highs the HiGHS instance to configure
inline void configure_highs(Highs& highs) {
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("log_to_console", false);
    highs.setOptionValue("solver", "simplex");
    highs.setOptionValue("simplex_strategy", 4);
    highs.setOptionValue("primal_feasibility_tolerance", 1e-7);
    highs.setOptionValue("dual_feasibility_tolerance", 1e-7);
    highs.setOptionValue("time_limit", 200);
}

// Loads the polytope A_eq x = b_eq, A_in x <= b_in, lb <= x <= ub into a HiGHS instance.
//
// The bounds lb, ub become column bounds, A_eq and A_in become rows.
// @tparam Point the point type
// @param P the polytope
// @param highs the HiGHS instance to load into, must be empty
template <typename Point>
void build_lp_model(ConstrainedPolytope<Point> const& P,
                    Highs& highs)
{   
    using MT = typename ConstrainedPolytope<Point>::MT;

    auto const& A_eq = P.getEqualities();
    auto const& b_eq = P.getEqualityRHS();
    auto const& A_in = P.getInequalities();
    auto const& b_in = P.getInequalityRHS();
    auto const& lb = P.getLowerBounds();
    auto const& ub = P.getUpperBounds();
    unsigned const d = P.getDimension();
    unsigned const m_eq = P.getNumEqualities();
    unsigned const m_in = P.getNumInequalities();

    // Adds the variables with their bounds
    for (unsigned j = 0; j < d; ++j) {
        double const l = (double)lb(j);
        double const u = (double)ub(j);
        highs.addVar(std::isfinite(l) ? l : -kHighsInf,
                     std::isfinite(u) ? u : kHighsInf);
    }

    // Adds a row of A with bounds lo <= a x <= hi
    auto add_row = [&](MT const& A, unsigned i, double lo, double hi) {
        std::vector<HighsInt> ids;
        std::vector<double> vals;

        for (typename MT::InnerIterator it(A, i); it; ++it) {
            ids.push_back((HighsInt)it.col());
            vals.push_back((double)it.value());
        }

        highs.addRow(lo, hi, (HighsInt)ids.size(), ids.data(), vals.data());
    };

    // Adds the equalities
    for (unsigned i = 0; i < m_eq; ++i)
        add_row(A_eq, i, (double)b_eq(i), (double)b_eq(i));

    // Adds the inequalities
    for (unsigned i = 0; i < m_in; ++i)
        add_row(A_in, i, -kHighsInf, (double)b_in(i));
}

// Extracts the current model of a HiGHS instance as a ConstrainedPolytope.
//
// Rows with equal lower and upper bounds become rows of A_eq. Every other row
// lo <= a x <= hi becomes up to two rows of A_in: a x <= hi and -a x <= -lo, one
// per finite side. Column bounds become lb, ub.
// @tparam Point the point type
// @param highs the HiGHS instance to read from
// @param tol the numerical tolerance
// @return the polytope extracted from the model
template <typename Point>
ConstrainedPolytope<Point> extract_polytope(Highs const& highs,
                                            double tol = 1e-7)
{
    using NT = typename ConstrainedPolytope<Point>::NT;
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;
    using Triplet = typename ConstrainedPolytope<Point>::Triplet;

    HighsLp const& lp = highs.getLp();
    unsigned const d = (unsigned)lp.num_col_;
    unsigned const num_rows = (unsigned)lp.num_row_;

    HighsSparseMatrix A = lp.a_matrix_;
    A.ensureRowwise();
    
    std::vector<Triplet> eq_triplets, in_triplets;
    std::vector<NT> eq_rhs, in_rhs;

    // Appends row i, scaled by sign to a system
    auto push_row = [&](std::vector<Triplet>& trip, std::vector<NT>& rhs, HighsInt i, double sign, double b) {
        int const row = (int)rhs.size();
        for (HighsInt p = A.start_[i]; p < A.start_[i+1]; ++p)
            trip.emplace_back(row, (int)A.index_[p], NT(sign*A.value_[p]));
        rhs.push_back(NT(sign*b));
    };

    for (unsigned i = 0; i < num_rows; ++i) {
        double const lo = lp.row_lower_[i];
        double const hi = lp.row_upper_[i];
        
        if (std::isfinite(lo) && std::isfinite(hi) && std::abs(hi-lo) < tol) {
            push_row(eq_triplets, eq_rhs, i, 1.0, hi);
            continue;
        }

        if (std::isfinite(hi)) push_row(in_triplets, in_rhs, i, 1.0, hi);
        if (std::isfinite(lo)) push_row(in_triplets, in_rhs, i, -1.0, lo);
    }

    unsigned const m_eq = (unsigned)eq_rhs.size();
    unsigned const m_in = (unsigned)in_rhs.size();

    MT A_eq(m_eq, d), A_in(m_in, d);
    A_eq.setFromTriplets(eq_triplets.begin(), eq_triplets.end());
    A_in.setFromTriplets(in_triplets.begin(), in_triplets.end());

    VT b_eq = Eigen::Map<VT const>(eq_rhs.data(), m_eq);
    VT b_in = Eigen::Map<VT const>(in_rhs.data(), m_in);

    NT const INF = std::numeric_limits<NT>::infinity();
    VT lb(d), ub(d);
    for (unsigned j = 0; j < d; ++j) {
        lb(j) = std::isfinite(lp.col_lower_[j]) ? NT(lp.col_lower_[j]) : -INF;
        ub(j) = std::isfinite(lp.col_upper_[j]) ? NT(lp.col_upper_[j]) : INF;
    }

    return ConstrainedPolytope<Point>(d, A_eq, b_eq, A_in, b_in, lb, ub);
}
#endif