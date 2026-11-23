// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef REDUCE_POLYTOPE_HPP
#define REDUCE_POLYTOPE_HPP

#include <cmath>
#include <vector>
#include <limits>
#include "Highs.h"
#include "convex_bodies/constrained_polytope.hpp"
#include "preprocess/constrained_polytope/highs_polytope.hpp"

// Configuration for reducing a polytope.
struct ReduceConfig {
    // A slack at or below this counts as zero.
    double slack_tol = 1e-7;
};

// Result of reducing a polytope
template <typename Point>
struct ReduceResult {
    // The reduced polytope.
    ConstrainedPolytope<Point> polytope;

    // The number of A_in rows found to be implicit equalities.
    unsigned implicit_rows = 0;
    
    // The number of variables fixed.
    unsigned fixed_vars = 0;

    // The number of LPs that were retried.
    unsigned retried_lps = 0;

    // The number of LPs solved.
    unsigned solved_lps = 0;

    // Tracks if the reduction algorithm ran.
    bool valid = false;
};

// Reduces a ConstrainedPolytope to its affine hull by finding the inequality rows
// and bounds that hold with equality at every feasible point, by moving them into A_eq.
//
// Every candidate r (a row of A_in or a finite bound) gets a slack s_r in [0, 1] scaled
// by the norm of the row:
//
//      max sum s_r, 
//      a_i x + ||a_i|| s_r <= b_i  rows of A_in
//      x_j+s_r <= ub(j)            upper bounds
//      -x_j+s_r <= -lb(j)          lower bounds
//
// A candidate with positive slack can be strict, its slack is fixed at 0 and it leaves
// the objective. The remaining ones are exactly the implicit equalities.
//
// The result is invalid only if the polytope is infeasible. If an LP fails even after
// a retry, the input polytope is returned unchanged, which is still correct.
// @tparam Point the point type
// @param P the input polytope
// @param config the config file
// @return the result
template <typename Point>
ReduceResult<Point> reduce_polytope(ConstrainedPolytope<Point> const& P,
                                    ReduceConfig const& config = ReduceConfig{})
{
    using NT = typename ConstrainedPolytope<Point>::NT;
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;
    using Triplet = typename ConstrainedPolytope<Point>::Triplet;

    ReduceResult<Point> res;

    MT const& A_eq = P.getEqualities();
    VT const& b_eq = P.getEqualityRHS();
    MT const& A_in = P.getInequalities();
    VT const& b_in = P.getInequalityRHS();
    VT const& lb = P.getLowerBounds();
    VT const& ub = P.getUpperBounds();
    unsigned const d = P.getDimension();
    unsigned const m_eq = P.getNumEqualities();
    unsigned const m_in = P.getNumInequalities();

    // Groups together bounds
    enum Kind { ROW, UPPER, LOWER };

    // A candidate is a constraint that may secretly be an equality
    struct Candidate { Kind kind; unsigned index; };

    std::vector<Candidate> cands;
    for (unsigned i = 0; i < m_in; ++i) cands.push_back({ROW, i});
    for (unsigned j = 0; j < d; ++j) {
        if (std::isfinite((double)ub(j))) cands.push_back({UPPER, j});
        if (std::isfinite((double)lb(j))) cands.push_back({LOWER, j});
    }

    unsigned const m = (unsigned)cands.size();

    // Nothing to reduce
    if (m == 0) {
        res.polytope = P;
        res.valid = true;
        return res;
    }

    Highs highs;
    configure_highs(highs);

    // Normals vars
    for (unsigned j = 0; j < d; ++j)
        highs.addVar(-kHighsInf, kHighsInf);

    // Adds the slack vars
    for (unsigned r = 0; r < m; ++r) {
        highs.addVar(0.0, 1.0);
        highs.changeColCost((HighsInt)(d+r), 1.0);
    }

    // Adds the equalities
    for (unsigned i = 0; i < m_eq; ++i) {
        std::vector<HighsInt> ids;
        std::vector<double> vals;
        for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
            ids.push_back((HighsInt)it.col());
            vals.push_back((double)it.value());
        }
        highs.addRow((double)b_eq(i), (double)b_eq(i), (HighsInt)ids.size(), ids.data(), vals.data());
    }

    // Adds one row per candidate
    for (unsigned r = 0; r < m; ++r) {
        HighsInt const s = (HighsInt)(d+r);
        Candidate const c = cands[r];

        if (c.kind == ROW) {
            // a_i x + ||a_i|| s_r <= b_i
            std::vector<HighsInt> ids;
            std::vector<double> vals;

            for (typename MT::InnerIterator it(A_in, c.index); it; ++it) {
                ids.push_back((HighsInt)it.col());
                vals.push_back((double)it.value());
            }
            ids.push_back(s);
            vals.push_back((double)A_in.row(c.index).norm());
            highs.addRow(-kHighsInf, (double)b_in(c.index), (HighsInt)ids.size(), ids.data(), vals.data());
        } else {
            HighsInt ids[2] = {(HighsInt)c.index, s};
            double const sign = (c.kind == UPPER) ? 1.0 : -1.0;
            double vals[2] = {sign, 1.0};
            double const rhs = (c.kind == UPPER) ? (double)ub(c.index) : -(double)lb(c.index);
            highs.addRow(-kHighsInf, rhs, 2, ids, vals);
        }
    }

    highs.changeObjectiveSense(ObjSense::kMaximize);

    std::vector<bool> active(m, true);
    unsigned n_active = m;
    while (n_active > 0) {
        highs.run();
        ++res.solved_lps;

        // The slacks are at most 1 and x is not in the objective, so the LP is never
        // unbounded: a status other than optimal or infeasible is a solver failure
        auto infeasible = [](HighsModelStatus s) {
            return s == HighsModelStatus::kInfeasible || s == HighsModelStatus::kUnboundedOrInfeasible;
        };

        auto s = highs.getModelStatus();
        if (s != HighsModelStatus::kOptimal && !infeasible(s)) {
            ++res.retried_lps;
            highs.clearSolver();
            highs.run();
            ++res.solved_lps;
            s = highs.getModelStatus();
        }

        if (infeasible(s)) return res;

        if (s != HighsModelStatus::kOptimal) {
            // Keeps the input, the active candidates are not proven implicit yet
            res.polytope = P;
            res.implicit_rows = 0;
            res.fixed_vars = 0;
            res.valid = true;
            return res;
        }

        auto const& sol = highs.getSolution().col_value;
        unsigned released = 0;

        for (unsigned r = 0; r < m; ++r) {
            if (!active[r] || sol[d+r] <= config.slack_tol) continue;

            active[r] = false;
            --n_active;
            ++released;

            highs.changeColBounds((HighsInt)(d+r), 0.0, 0.0);
            highs.changeColCost((HighsInt)(d+r), 0.0);
        }

        if (released == 0) break;
    }
    
    // Builds the output
    NT const INF = std::numeric_limits<NT>::infinity();
    std::vector<Triplet> eq_trip, in_trip;
    std::vector<NT> eq_rhs, in_rhs;
    VT lb_new = lb, ub_new = ub;
    std::vector<bool> fixed(d, false);

    for (unsigned i = 0; i < m_eq; ++i) {
        for (typename MT::InnerIterator it(A_eq, i); it; ++it)
            eq_trip.emplace_back((int)eq_rhs.size(), (int)it.col(), it.value());
        eq_rhs.push_back(b_eq(i));
    }

    for (unsigned r = 0; r < m; ++r) {
        Candidate const c = cands[r];

        if (c.kind == ROW) {
            bool const implicit = active[r];
            auto& trip = implicit ? eq_trip : in_trip;
            auto& rhs = implicit ? eq_rhs : in_rhs;

            for (typename MT::InnerIterator it(A_in, c.index); it; ++it)
                trip.emplace_back((int)rhs.size(), (int)it.col(), it.value());
            rhs.push_back(b_in(c.index));

            if (implicit) ++res.implicit_rows;
            continue;
        }

        if (!active[r] || fixed[c.index]) continue;

        unsigned const j = c.index;
        NT const val = (c.kind == UPPER) ? ub(j) : lb(j);

        eq_trip.emplace_back((int)eq_rhs.size(), (int)j, NT(1));
        eq_rhs.push_back(val);

        lb_new(j) = -INF;
        ub_new(j) = INF;
        fixed[j] = true;
        ++res.fixed_vars;
    }

    MT A_eq_new((Eigen::Index)eq_rhs.size(), d), A_in_new((Eigen::Index)in_rhs.size(), d);
    A_eq_new.setFromTriplets(eq_trip.begin(), eq_trip.end());
    A_in_new.setFromTriplets(in_trip.begin(), in_trip.end());

    VT b_eq_new = Eigen::Map<VT const>(eq_rhs.data(), (Eigen::Index)eq_rhs.size());
    VT b_in_new = Eigen::Map<VT const>(in_rhs.data(), (Eigen::Index)in_rhs.size());
    res.polytope = ConstrainedPolytope<Point>(d, A_eq_new, b_eq_new, A_in_new, b_in_new, lb_new, ub_new);
    res.valid = true;
    return res;
}
#endif