// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONSTRAINED_POLYTOPE_SIMPLIFIER_HPP
#define CONSTRAINED_POLYTOPE_SIMPLIFIER_HPP

#include <ctime>
#include <algorithm>
#include <string>
#include "convex_bodies/constrained_polytope.hpp"
#include "preprocess/constrained_polytope/reduce_polytope.hpp"
#include "preprocess/constrained_polytope/interior_point.hpp"
#include "preprocess/constrained_polytope/clarkson.hpp"
#include "preprocess/constrained_polytope/presolve.hpp"

// Configuration for the simplification pipeline.
struct SimplifierConfig {
    // Configuration for dimension fixing.
    ReduceConfig reduce;

    // Configuration for dependent row removal.
    DependentRowsRemovalConfig dependent_rows;

    // Configuration for Clarkson's redundancy removal.
    ClarksonConfig clarkson;

    // The minimum slack for the interior point to count.
    double interior_tol = 1e-9;

    // Tracks if the reducer should be executed.
    bool fix_dimensions = true;

    // Tracks if presolve should happen.
    bool presolve = true;

    // Tracks if dependent rows should be removed.
    bool remove_dependent_rows = true;

    // The maximum number of polytope reduction passes. When
    // no interior point is found after a pass, dimension fixing
    // is attempted again.
    unsigned max_reduction_passes = 2;
};

// The outcome of a simplification run.
enum class SimplifierStatus {
    OK,
    PRESOLVE_FAILED,
    REMOVE_DEPENDENT_ROWS_FAILED,
    INFEASIBLE,
    NO_INTERIOR_POINT,
    CLARKSON_FAILED
};

// Returns the name of a simplifier status.
inline std::string status_to_string(SimplifierStatus s) {
    switch (s) {
        case SimplifierStatus::OK: return "OK";
        case SimplifierStatus::PRESOLVE_FAILED: return "PRESOLVE_FAILED";
        case SimplifierStatus::REMOVE_DEPENDENT_ROWS_FAILED: return "REMOVE_DEPENDENT_ROWS_FAILED";
        case SimplifierStatus::INFEASIBLE: return "INFEASIBLE";
        case SimplifierStatus::NO_INTERIOR_POINT: return "NO_INTERIOR_POINT";
        case SimplifierStatus::CLARKSON_FAILED: return "CLARKSON_FAILED";
    }
    return "UNKNOWN";
}

// The cpu time of every stage in seconds.
struct SimplifierTimes {
    double presolve = 0.0;
    double reduce = 0.0;
    double dependent_rows = 0.0;
    double interior_point = 0.0;
    double clarkson = 0.0;

    // @return the CPU time of the whole run
    double total() const {
        return presolve+reduce+dependent_rows+interior_point+clarkson;
    }
};

// Statistics about presolve.
struct PresolveStats {
    // Number of variables pinned by A_eq x = b_eq.
    unsigned pinned_vars = 0;

    // The rank of A_eq.
    unsigned rank = 0;
};

// Statistics about the reduce stage.
struct ReduceStats {
    // The number of passes.
    unsigned passes = 0;

    // The number of variables fixed.
    unsigned fixed_vars = 0;

    // The number of LPs solved.
    unsigned solved_lps = 0;
};

// Statistics about the dependent row removal stage.
struct DependentRowsStats {
    // The rank of A_eq.
    unsigned rank = 0;

    // The number of rows removed.
    unsigned removed_rows = 0;

    // The largest violation of a dropped row observed.
    double max_residual = 0.0;
};

// Statistics about the Clarkson stage.
struct ClarksonStats {
    // The number of LPs solved.
    unsigned solved_lps = 0;

    // The number of LPs retried.
    unsigned retried_lps = 0;

    // The number of rows kept untested after repeated failures.
    unsigned untested_rows = 0;
};

// Result of the simplification pipeline.
template <typename Point> 
struct SimplifierResult {
   // The simplified polytope.
   ConstrainedPolytope<Point> polytope; 
   
   // Number of LPs solved.
   unsigned solved_lps = 0;

   // The outcome of simplification
   SimplifierStatus status = SimplifierStatus::OK;

   // The CPU time of every stage in seconds.
   SimplifierTimes times;

   // Statistics about presolve.
   PresolveStats presolve;

   // Statistics about reduce.
   ReduceStats reduce;

   // The interior point found.
   InteriorPoint<typename ConstrainedPolytope<Point>::VT> interior_point;

   // Statistics about dependent row removal.
   DependentRowsStats dependent_rows;

   // Statistics about Clarkson.
   ClarksonStats clarkson;
};

// Simplifies a ConstrainedPolytope by running reducing it, finding an interior point
// and removing redundant inequality rows using Clarkson's algorithm.
//
// Pipeline:
//      1. presolve_pinned_vars, frees bounds of variables pinned by A_eq x = b_eq
//      2. reduce_polytope, fixes degenerate dimensions
//      3. remove_dependent_rows, removes dependent rows
//      4. find_interior_point, finds a strictly interior point
//      5. redundancy_removal_clarkson, removes redundant rows from A_in
// @tparam Point the point type
// @param P the input polytope
// @param config the simplification configuration
// @return the result
template <typename Point>
SimplifierResult<Point> simplify(ConstrainedPolytope<Point> const& P,
                                 SimplifierConfig const& config = SimplifierConfig{})
{
    SimplifierResult<Point> res;
    ConstrainedPolytope<Point> Ps = P;

    auto fail = [&](SimplifierStatus s) {
        res.status = s;
        res.polytope = Ps;
        return res;
    };

    auto elapsed = [](std::clock_t c0) {
        return double(std::clock()-c0)/CLOCKS_PER_SEC;
    };

    // Stage 1: relax pinned variables by A_eq = b_eq (presolve)
    if (config.presolve) {
        std::clock_t c0 = std::clock();
        auto pres = presolve_pinned_vars(Ps);
        res.times.presolve = elapsed(c0);

        res.presolve.pinned_vars = (unsigned)pres.pinned.size();
        res.presolve.rank = pres.rank;

        if (!pres.valid) return fail(SimplifierStatus::PRESOLVE_FAILED);
        Ps = pres.polytope;
    }

    for (unsigned pass = 1; pass <= config.max_reduction_passes; ++pass) {
        // Stage 2: fix degenerate dimensions
        if (config.fix_dimensions) {
            std::clock_t c0 = std::clock();
            auto rres = reduce_polytope(Ps, config.reduce);
            res.times.reduce += elapsed(c0);

            res.reduce.passes += 1;
            res.reduce.fixed_vars += rres.fixed_vars;
            res.reduce.solved_lps += rres.solved_lps;
            res.solved_lps += rres.solved_lps;

            if (!rres.valid) return fail(SimplifierStatus::INFEASIBLE);
            Ps = rres.polytope;
        }

        // Stage 3: dependent equality rows
        if (config.remove_dependent_rows) {
            std::clock_t c0 = std::clock();
            auto dres = remove_dependent_rows(Ps, config.dependent_rows);
            res.times.dependent_rows += elapsed(c0);

            res.dependent_rows.rank = dres.rank;
            res.dependent_rows.removed_rows += dres.removed_rows;
            res.dependent_rows.max_residual = std::max(res.dependent_rows.max_residual, dres.max_residual);

            if (!dres.valid) return fail(SimplifierStatus::REMOVE_DEPENDENT_ROWS_FAILED);
            Ps = dres.polytope;
        }

        // Stage 4: Interior point
        std::clock_t c0 = std::clock();
        res.interior_point = find_interior_point(Ps, config.interior_tol);
        res.times.interior_point += elapsed(c0);
        res.solved_lps += 1;

        if (!res.interior_point.valid) {
            if (config.fix_dimensions && pass < config.max_reduction_passes) {
                continue; // try reducing the polytope again
            }
            return fail(SimplifierStatus::NO_INTERIOR_POINT);
        }

        break; // interior point found
    }

    // Stage 5: Clarkson redundancy removal
    std::clock_t c0 = std::clock();
    auto cres = redundancy_removal_clarkson(Ps, res.interior_point.point, config.clarkson);
    res.times.clarkson = elapsed(c0);

    res.clarkson.solved_lps += cres.solved_lps;
    res.clarkson.retried_lps += cres.retried_lps;
    res.clarkson.untested_rows += cres.untested_rows;
    res.solved_lps += cres.solved_lps;

    if (!cres.valid) return fail(SimplifierStatus::CLARKSON_FAILED);

    res.polytope = cres.polytope;
    return res;
}

#endif