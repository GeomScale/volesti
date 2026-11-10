// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CLARKSON_HPP
#define CLARKSON_HPP

#include <random>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <utility>
#include "Highs.h"
#include "convex_bodies/constrained_polytope.hpp"
#include "preprocess/constrained_polytope/highs_polytope.hpp"
#include "preprocess/constrained_polytope/ray_shoot.hpp"

// Configuration for Clarkson's redundancy removal.
struct ClarksonConfig {
    // The gap by which Clarkson relaxes a bound.
    double relaxation_gap = 1.0;

    // The tolerance for declaring a row redundant.
    double facet_tol =1e-7;

    // Ray directions below this are treated as zero.
    double ray_tol = 1e-9;

    // The number of times an LP for one row may fail before it is kept.
    unsigned failed_iter_count = 2;

    // The seed for the random row picker.
    unsigned seed = 1;
};

// Stores the result of Clarkson's algorithm.
template <typename Point>
struct ClarksonResult {
    // The simplified polytope with redundant bounds removed.
    ConstrainedPolytope<Point> polytope;

    // The number of LPs solved.
    unsigned solved_lps = 0;

    // The number of LPs retried.
    unsigned retried_lps = 0;

    // The number of rows kept untested after repeated failures.
    unsigned untested_rows = 0;

    // True if the algorithm ran successfully.
    bool valid = false;
};

// Tests whether row k of A_in is redundant given the current essential set already
// enforced in the HiGHS model.
//
// The row is temporarily relaxed by the relaxation_gap and the LP maximizes a_k^T x.
// @tparam Point the point type
// @param highs the current LP model with the current essential constraints enforced
// @param A_in the inequality matrix
// @param b_in the inequality rhs
// @param k the row to test
// @param m_eq number of equalities
// @param m_in number of inequalities
// @param d the number of variables
// @param config the Clarkson configuration
// @param solved_lps counter incremented for every HiGHS run
// @param retried_lps counter incremented for every failed HiGHS run
// @param solved set to true if the LP solved to optimality
// @return whether the row is redundant and the LP optimum x*
template <typename Point>
std::pair<typename ConstrainedPolytope<Point>::VT, bool>
test_redundancy(Highs& highs, 
                typename ConstrainedPolytope<Point>::MT const& A,
                typename ConstrainedPolytope<Point>::VT const& b,
                unsigned k,
                unsigned m_eq,
                unsigned m_in,
                std::vector<std::pair<unsigned, bool>> const& bounds,
                unsigned d,
                ClarksonConfig const& config,
                unsigned& solved_lps,
                unsigned& retried_lps,
                bool& solved)
{
    using NT = typename ConstrainedPolytope<Point>::NT;
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;

    double const relaxed = (double)b(k)+config.relaxation_gap;
    bool const is_bound = (k >= m_in);
    HighsInt test_row = -1;

    HighsInt index = 0;
    double old_lo = 0.0;
    double old_hi = 0.0;

    if (is_bound) {
        auto [j, is_upper] = bounds[k-m_in];
        index = (HighsInt)j;
        old_lo = highs.getLp().col_lower_[index];
        old_hi = highs.getLp().col_upper_[index];

        if (is_upper) {
            highs.changeColBounds(index, old_lo, relaxed);
        } else {
            highs.changeColBounds(index, -relaxed, old_hi);
        }
    } else {
        index = (HighsInt)(m_eq+k);
        old_lo = highs.getLp().row_lower_[index];
        old_hi = highs.getLp().row_upper_[index];
        highs.changeRowBounds(index, -kHighsInf, relaxed);
    }

    // Maximizes the constraint
    for (typename MT::InnerIterator it(A, k); it; ++it)
        highs.changeColCost((HighsInt)it.col(), (double)it.value());
    
    highs.changeObjectiveSense(ObjSense::kMaximize);

    highs.run();
    ++solved_lps;

    if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
        highs.clearSolver();
        highs.run();
        ++solved_lps;
        ++retried_lps;
    }

    solved = (highs.getModelStatus() == HighsModelStatus::kOptimal);

    VT x_star(d);
    bool redundant = false;

    if (solved) {
        auto const& sol = highs.getSolution().col_value;
        for (unsigned j = 0; j < d; ++j)
            x_star(j) = (NT)sol[j];

        double ax = (double)A.row(k).dot(x_star);
        double const obj = highs.getObjectiveValue();
        redundant = (std::max(ax, obj)) <= (double)b(k)+config.facet_tol;
    }

    // Restores the model
    if (is_bound) {
        highs.changeColBounds(index, old_lo, old_hi);
    } else {
        highs.changeRowBounds(index, old_lo, old_hi);
    }

    for (typename MT::InnerIterator it(A, k); it; ++it)
        highs.changeColCost((HighsInt)it.col(), 0.0);

    return {x_star, redundant};
}


// Builds the base LP for Clarkson's algorithm:
//
//      A_eq x = b_eq
//
// @tparam Point the point type
// @param P the polytope
// @param highs an empty HiGHS instance to load into
template <typename Point>
void build_clarkson_lp(ConstrainedPolytope<Point> const& P, Highs& highs)
{
    using MT = typename ConstrainedPolytope<Point>::MT;

    auto const& A_eq = P.getEqualities();
    auto const& b_eq = P.getEqualityRHS();
    unsigned const d = P.getDimension();
    auto const& A_in = P.getInequalities();
    unsigned const m_in = P.getNumInequalities();
    unsigned const m_eq = P.getNumEqualities();

    for (unsigned j = 0; j < d; ++j)
        highs.addVar(-kHighsInf, kHighsInf);

    for (unsigned i = 0; i < m_eq; ++i) {
        std::vector<HighsInt> ids;
        std::vector<double> vals;
        for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
            ids.push_back((HighsInt)it.col());
            vals.push_back((double)it.value());
        }

        highs.addRow((double)b_eq(i), (double)b_eq(i), (HighsInt)ids.size(),
                     ids.data(), vals.data());
    }

    for (unsigned i = 0; i < m_in; ++i) {
        std::vector<HighsInt> ids;
        std::vector<double> vals;
        for (typename MT::InnerIterator it(A_in, i); it; ++it) {
            ids.push_back((HighsInt)it.col());
            vals.push_back((double)it.value());
        }

        highs.addRow(-kHighsInf, kHighsInf, (HighsInt)ids.size(),
                     ids.data(), vals.data());
    }
}

// Removes redundant rows from the inequality system A_in x <= b_in using Clarkson's
// algorithm.
//
// The algorithm starts with only the equality constraints in the LP and adds inequality rows
// one at a time as they are proved essential.
//
// @tparam Point the point type
// @param P the polytope whose inequalities are to be simplified
// @param z a strictly interior point of P
// @param config the algorithm configuration
// @return the result
template <typename Point>
ClarksonResult<Point> redundancy_removal_clarkson(ConstrainedPolytope<Point> const& P,
                                                  typename ConstrainedPolytope<Point>::VT const& z,
                                                  ClarksonConfig const& config = ClarksonConfig{})
{
    using NT = typename ConstrainedPolytope<Point>::NT;
    using VT = typename ConstrainedPolytope<Point>::VT;
    using MT = typename ConstrainedPolytope<Point>::MT;
    using Triplet = typename ConstrainedPolytope<Point>::Triplet;

    auto const& A_eq = P.getEqualities();
    auto const& b_eq = P.getEqualityRHS();
    auto const& A_in = P.getInequalities();
    auto const& b_in = P.getInequalityRHS();
    auto const& lb = P.getLowerBounds();
    auto const& ub = P.getUpperBounds();
    unsigned const d = P.getDimension();
    unsigned const m_eq = P.getNumEqualities();
    unsigned const m_in = P.getNumInequalities();

    // Stacks the finite bounds under A_in
    std::vector<Triplet> triplets;
    std::vector<NT> rhs;
    std::vector<std::pair<unsigned, bool>> bounds;

    for (unsigned i = 0; i < m_in; ++i) {
        for (typename MT::InnerIterator it(A_in, i); it; ++it)
            triplets.emplace_back((int)i, (int)it.col(), it.value());
        rhs.push_back(b_in(i));
    }

    for (unsigned j = 0; j < d; ++j) {
        if (std::isfinite((double)ub(j))) {
            triplets.emplace_back((int)rhs.size(), (int)j, NT(1));
            rhs.push_back(ub(j));
            bounds.emplace_back(j, true);
        } 
        if (std::isfinite((double)lb(j))) {
            triplets.emplace_back((int)rhs.size(), (int)j, NT(-1));
            rhs.push_back(-lb(j));
            bounds.emplace_back(j, false);
        } 
    }

    unsigned const m = (unsigned)rhs.size();

    MT A(m, d);
    A.setFromTriplets(triplets.begin(), triplets.end());
    VT b = Eigen::Map<VT const>(rhs.data(), m);

    // Same polytope, but bounds are in A_in, b_in (useful for rayshooting)
    ConstrainedPolytope<Point> Pc(d, A_eq, b_eq, A, b);

    ClarksonResult<Point> res;

    Highs highs;
    configure_highs(highs);
    highs.setOptionValue("presolve", "off");
    build_clarkson_lp(P, highs);

    // Enforces row i by adding it permanetly to the HiGHS model
    auto enforce = [&](unsigned i) {
        if (i >= m_in) {
            auto [j, is_upper] = bounds[i-m_in];
            double lo = highs.getLp().col_lower_[j];
            double hi = highs.getLp().col_upper_[j];

            if (is_upper) {
                hi = (double)b(i);
            } else {
                lo = -(double)b(i);
            }

            highs.changeColBounds((HighsInt)j, lo, hi);
            return;
        }

        highs.changeRowBounds((HighsInt)(m_eq+i), -kHighsInf, (double)b_in(i));
    };

    // The candidate rows whose redundancy is unknown
    std::vector<unsigned> J;

    // The position of row i in J, -1 if not in J
    std::vector<int> pos(m, -1);

    auto J_insert = [&](unsigned i) {
        pos[i] = (int)J.size();
        J.push_back(i);
    };

    auto J_erase = [&](unsigned i) -> bool {
        int p = pos[i];
        if (p < 0) return false;
        unsigned last = J.back();
        J[p] = last;
        pos[last] = p;
        J.pop_back();
        pos[i] = -1;
        return true;
    };

    for (unsigned i = 0; i < m; ++i)
        J_insert(i);

    std::vector<unsigned> I;
    std::mt19937 rng(config.seed);
    std::vector<unsigned> fail_count(m, 0);

    while (!J.empty()) {
        std::uniform_int_distribution<std::size_t> pick(0, J.size()-1);
        unsigned k = J[pick(rng)];

        bool solved = false;
        auto [x_star, redundant] = test_redundancy<Point>(
            highs, A, b, k, m_eq, m_in, bounds, d, config, res.solved_lps, res.retried_lps,
            solved);

        if (!solved) {
            if (++fail_count[k] > config.failed_iter_count) {
                res.untested_rows = (unsigned)J.size();
                for (unsigned i : J) I.push_back(i);
                break;
            }
            continue;
        }

        if (redundant) {
            J_erase(k);
            continue;
        }

        RayShoot rs = ray_shoot(Pc, J, z, x_star-z, config.ray_tol);

        if (!rs.hit || !J_erase(rs.row)) {
            I.push_back(k);
            enforce(k);
            J_erase(k);
        } else {
            I.push_back(rs.row);
            enforce(rs.row);
        }
    }

    std::vector<bool> keep(m, false);
    for (unsigned i : I) keep[i] = true;

    std::vector<Triplet> triplets_out;
    std::vector<NT> rhs_out;

    for (unsigned i = 0; i < m_in; ++i) {
        if (!keep[i]) continue;
        for (typename MT::InnerIterator it(A_in, i); it; ++it) 
            triplets_out.emplace_back((int)rhs_out.size(), (int)it.col(), it.value());
        rhs_out.push_back(b_in(i));
    }

    MT A_in_new((Eigen::Index)rhs_out.size(), d);
    A_in_new.setFromTriplets(triplets_out.begin(), triplets_out.end());
    VT b_in_new = Eigen::Map<VT const>(rhs_out.data(), (Eigen::Index)rhs_out.size());

    NT const INF = std::numeric_limits<NT>::infinity();
    VT lb_new = VT::Constant(d, -INF);
    VT ub_new = VT::Constant(d, INF);

    for (unsigned i = m_in; i < m; ++i) {
        if (!keep[i]) continue;
        auto [j, is_upper] = bounds[i-m_in];
        if (is_upper) {
            ub_new(j) = ub(j);
        } else {
            lb_new(j) = lb(j);
        }
    }

    res.polytope = ConstrainedPolytope<Point>(d, A_eq, b_eq, A_in_new, b_in_new, lb_new, ub_new);
    res.valid = true;
    return res;
}
#endif