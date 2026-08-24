// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_POLYORACLES_HPP
#define METABOLIC_POLYORACLES_HPP

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include "Highs.h"
#include "convex_bodies/metabolic_polytope.hpp"
#include "lp_oracles/lp_oracle_options.hpp"

// Checks whether P is contained in Q, by applying every constraint of Q on P.
// @tparam Point the point type of the polytope
// @param P the polytope that should be contained
// @param Q the polytope whose constraints are tested
// @param tol the constraint violation tolerance
// @param opts the LP solver options
// @return whether P is contained in Q, and whether every LP was solved
template <typename Point>
LPOracleResult<bool> is_contained_in(MetabolicPolytope<Point> const& P,
                                     MetabolicPolytope<Point> const& Q,
                                     double tol = 1e-3,
                                     LPOracleOptions const& opts = nullptr)
{   
    typedef typename MetabolicPolytope<Point>::MT MT;
    typedef typename MetabolicPolytope<Point>::VT VT;

    LPOracleResult<bool> res;
    res.solved = true;

    if (P.getDimension() != Q.getDimension()) return res;

    auto polytope_optimum = [&](Highs & highs, bool maximize, double & value) {
        highs.changeObjectiveSense(maximize ? ObjSense::kMaximize
                                            : ObjSense::kMinimize);
        
        highs.run();
        auto st = highs.getModelStatus();

        if (st == HighsModelStatus::kOptimal) {
            value = highs.getObjectiveValue();
            return true;
        }

        if (st == HighsModelStatus::kUnbounded) {
            value = maximize ? std::numeric_limits<double>::infinity()
                            : -std::numeric_limits<double>::infinity();
            
            return true;
        }

        #ifdef VOLESTI_DEBUG
        std::cout << "Could not solve the LP for metabolic polytope containment "
                  << ", highs returned status code " 
                  << (int)highs.getModelStatus()
                  << std::endl;
        #endif

        res.solved = false;
        return false;
    };

    unsigned d = P.getDimension();
    VT const& b_l = Q.getLowerBounds();
    VT const& b_u = Q.getUpperBounds();

    Highs highs;
    build_lp_model(P, highs);
    lp_oracles_configure_highs(highs, opts);

    for (unsigned k = 0; k < d; ++k) {
        bool has_up = !std::isinf((double)b_u(k));
        bool has_lo = !std::isinf((double)b_l(k));

        if (!has_up && !has_lo) continue;

        highs.changeColCost((HighsInt)k, 1.0);
        double v;

        if (has_up && (!polytope_optimum(highs, true, v) || v > (double)b_u(k)+tol))
            return res;

        if (has_lo && (!polytope_optimum(highs, false, v) || v < (double)b_l(k)-tol))
            return res;

        highs.changeColCost((HighsInt)k, 0.0);
    }

    MT const& A_eq = Q.getEqualities();
    VT const& b_eq = Q.getEqualityBounds();

    for (unsigned i = 0; i < (unsigned)A_eq.rows(); ++i) {
        std::vector<unsigned> cols;
        for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
            cols.push_back((unsigned)it.col());
            highs.changeColCost((HighsInt)it.col(), (double)it.value());
        }

        double rhs = (double)b_eq(i);
        double hi, lo;

        if (!polytope_optimum(highs, true, hi) || 
            !polytope_optimum(highs, false, lo) ||
            hi > rhs+tol || lo < rhs-tol) 
        {
            return res;
        }

        for (unsigned k : cols)
            highs.changeColCost((HighsInt)k, 0.0);
    }

    res.value = true;
    return res;
}

// Checks that two metabolic polytopes describe the same feasible region.
// @tparam Point the point type of the polytopes
// @param P the first polytope
// @param Q the second polytope
// @param tol the constraint violation tolerance
// @param opts the LP solver options
// @return whether P is equal to Q, and whether every LP was solved
template <typename Point>
LPOracleResult<bool> are_equal(MetabolicPolytope<Point> const& P,
                               MetabolicPolytope<Point> const& Q,
                               double tol = 1e-3,
                               LPOracleOptions const& opts = nullptr)
{   
    auto res_1 = is_contained_in(P, Q, tol, opts);
    auto res_2 = is_contained_in(Q, P, tol, opts);
    return {res_1.value && res_2.value, res_1.solved && res_2.solved};
}
#endif