// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONSTRAINED_POLYTOPE_LP_PARSER_HPP
#define CONSTRAINED_POLYTOPE_LP_PARSER_HPP

#include <stdexcept>
#include <string>
#include "Highs.h"
#include "convex_bodies/constrained_polytope.hpp"
#include "preprocess/constrained_polytope/highs_polytope.hpp"

// Parses the feasible region of a linear program into a ConstrainedPolytope:
//
//      A_eq x = b_eq (equality rows)
//      A_in x <= b_in (inequality rows)
//      lb <= x <= ub (variable bounds)
// The file is read by HiGHS, which chooses the format from the extension. The
// supported formats are MPS (.mps or .mps.gz) or CPLEX LP (.lp).
// @tparam Point the Point type
// @param model_path path to the model file
// @return the feasible region as a ConstrainedPolytope
template <typename Point>
ConstrainedPolytope<Point> parse_from_lp_model(std::string const& model_path) 
{
    Highs highs;
    highs.setOptionValue("output_flag", false);

    if (highs.readModel(model_path) == HighsStatus::kError)
        throw std::runtime_error("Cannot read the model "+model_path);

    // Throws an error if the model is ILP
    for (HighsVarType t : highs.getLp().integrality_)
        if (t != HighsVarType::kContinuous)
            throw std::runtime_error("The model "+model_path+" is not an LP");

    return extract_polytope<Point>(highs);
}
#endif