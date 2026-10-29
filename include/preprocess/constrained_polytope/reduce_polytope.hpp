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
#include "convex_bodies/constrained_polytope.hpp"
#include "preprocess/constrained_polytope/highs_polytope.hpp"

// Configuration for reducing a polytope.
struct ReduceConfig {
    // Variables whose observed range is below this are fixed.
    double reduce_tol = 1e-7;

    // Observed variation below this is treated as zero.
    double observe_tol = 1e-5;
};

// Result of reducing a polytope
template <typename Point>
struct ReduceResult {
    // The reduced polytope.
    ConstrainedPolytope<Point> polytope;

    // The number of variables fixed.
    unsigned fixed_vars = 0;

    // The number of LPs solved.
    unsigned solved_lps = 0;

    // Tracks if the reduction algorithm ran.
    bool valid = false;
};

// Fixes degenerate dimensions of a ConstrainedPolytope by identifying variables that take
// a single value over the feasible set and moving them into A_eq.
//
// The method builds a HiGHS LP from the polytope, runs one LP to collect observed
// variable values, and for each variable that shows no variation across the solutions
// seen so far it solves a max and min LP to confirm the range.
//
// The reduced polytope is returned with the fixed variables expressed as rows of A_eq.
// @tparam Point the point type
// @param P the input polytope
// @param config the config file
// @return the result
template <typename Point>
ReduceResult<Point> reduce_polytope(ConstrainedPolytope<Point> const& P,
                                    ReduceConfig const& config = ReduceConfig{})
{
    unsigned const d = P.getDimension();

    ReduceResult<Point> res;

    Highs highs;
    configure_highs(highs);
    build_lp_model(P, highs);

    auto run_lp = [&]() -> bool {
        highs.run();
        ++res.solved_lps;

        if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
            highs.clearSolver();
            highs.run();
            ++res.solved_lps;
        }

        return highs.getModelStatus() == HighsModelStatus::kOptimal;
    };

    // Adds x_k = val to A_eq and frees the bounds of x_k
    auto fix_dimension = [&](unsigned k, double val) {
        if (std::abs(val) < config.reduce_tol) val = 0.0;

        HighsInt id = (HighsInt)k;
        double coeff = 1.0;
        highs.addRow(val, val, 1, &id, &coeff);
        highs.changeColBounds(id, -kHighsInf, kHighsInf);
        ++res.fixed_vars;
    };

    // Checks that this region describes a non empty and bounded polytope
    if (!run_lp()) return res;

    // Observed min/max values for each variables across all LP solutions seen
    std::vector<double> obs_min(d, std::numeric_limits<double>::infinity());
    std::vector<double> obs_max(d, -std::numeric_limits<double>::infinity());

    auto observe = [&]() {
        auto const& sol = highs.getSolution().col_value;
        for (unsigned j = 0; j < d; ++j) {
            if (sol[j] < obs_min[j]) obs_min[j] = sol[j];
            if (sol[j] > obs_max[j]) obs_max[j] = sol[j];
        }
    };

    observe();

    auto observe_variation = [&](unsigned k) {
        return std::abs(obs_max[k]-obs_min[k]) > config.observe_tol;
    };

    for (unsigned k = 0; k < d; ++k) {
        double const l = highs.getLp().col_lower_[k];
        double const u = highs.getLp().col_upper_[k];

        // Prunes if the bounds are too close
        if (l > -kHighsInf && u < kHighsInf && std::abs(u-l) < config.reduce_tol) {
            fix_dimension(k, (l+u)/2.0);
            continue;
        }

        // Prunes based on past observation
        if (observe_variation(k)) continue;

        highs.changeColCost((HighsInt)k, 1.0);
        highs.changeObjectiveSense(ObjSense::kMaximize);

        if (!run_lp()) {
            highs.changeColCost((HighsInt)k, 0.0);
            continue;
        }

        double max_val = highs.getObjectiveValue();
        observe();

        if (observe_variation(k)) {
            highs.changeColCost((HighsInt)k, 0.0);
            continue;
        }

        highs.changeObjectiveSense(ObjSense::kMinimize);

        if (!run_lp()) {
            highs.changeColCost((HighsInt)k, 0.0);
            continue;
        }

        double min_val = highs.getObjectiveValue();
        observe();

        highs.changeColCost((HighsInt)k, 0.0);

        if (std::abs(max_val-min_val) < config.reduce_tol) {
            fix_dimension(k, (max_val+min_val)/2.0);
        }
    }
    
    res.polytope = extract_polytope<Point>(highs);
    res.valid = true;
    return res;
}
#endif