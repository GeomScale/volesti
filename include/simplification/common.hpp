// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SIMPLIFICATION_COMMON_HPP
#define SIMPLIFICATION_COMMON_HPP

#include "convex_bodies/metabolic_polytope.h"
#include "Highs.h"

namespace simplification {
    struct Config {
        double facet_tolerance = 1e-7;
        double dim_tolerance = 1e-7;
        double numerical_tolerance = 1e-7;
        bool verbose = false;
        bool fix_dimensions = false;
    };

    template <typename Point>
    struct Result {
        // The simplified polytope
        MetabolicPolytope<Point> P;

        // Simplification statistics
        unsigned bounds_relaxed = 0;  // bounds relaxed
        unsigned dims_fixed = 0;      // inequalities converted to equalities

        // Tracks if simplification was successful
        bool success = true;
    };
}
#endif