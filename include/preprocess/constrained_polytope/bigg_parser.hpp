// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BIGG_PARSER_HPP
#define BIGG_PARSER_HPP

#include <nlohmann/json.hpp>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <string>
#include <limits>
#include "convex_bodies/constrained_polytope.hpp"

// Parses a BiGG JSON model into a ConstrainedPolytope of the form:
//
//      A_eq x = b_eq  (steady-state condition S x = 0)
//      lb <= x <= ub  (flux bounds)
//
// Each reaction becomes a variable, each metabolite an equality row. Flux
// bounds are stored as bounds on x so A_in is empty.
// @tparam Point the Point type used by the metabolic network
// @param jsn a parsed nlohmann::json object holding the model
// @return the metabolic network as a ConstrainedPolytope instance
template <typename Point>
ConstrainedPolytope<Point> construct_from_json(nlohmann::json const& jsn) 
{
    using Polytope = ConstrainedPolytope<Point>;
    using MT = typename Polytope::MT;
    using VT = typename Polytope::VT;
    using NT = typename Polytope::NT;
    using Triplet = typename Polytope::Triplet;

    const NT INF = std::numeric_limits<NT>::infinity();

    if (!jsn.contains("metabolites") || !jsn.contains("reactions"))
        throw std::runtime_error("Not BiGG file: missing `metabolites` or `reactions` field");
        
    auto const& metabolites = jsn.at("metabolites");
    auto const& reactions = jsn.at("reactions");
    unsigned m = (unsigned)metabolites.size(); // metabolite count
    unsigned n = (unsigned)reactions.size();   // reaction count

    // Maps the metabolites to the integers in [0,m)
    std::unordered_map<std::string, unsigned> metabolite_index;
    {
        unsigned i = 0;
        for (auto const& metabolite : metabolites) {
            std::string met_id = metabolite.at("id").get<std::string>();
            metabolite_index.emplace(met_id, i++);
        }
    }

    // Builds the stoichiometric matrix A_eq and the box bounds lb, ub.
    std::vector<Triplet> eq_triplets;
    VT lb(n), ub(n);

    unsigned j = 0;

    for (auto const& reaction : reactions) {
        lb(j) = reaction.contains("lower_bound") ? reaction.at("lower_bound").get<NT>() : -INF;
        ub(j) = reaction.contains("upper_bound") ? reaction.at("upper_bound").get<NT>() : INF;

        // Stoichiometric coefficients.
        if (reaction.contains("metabolites")) {
            for (auto it = reaction.at("metabolites").begin(); it != reaction.at("metabolites").end(); ++it) {
                auto found = metabolite_index.find(it.key());
                if (found == metabolite_index.end()) // metabolite wasn't found, file is corrupt
                    throw std::runtime_error("Reaction references unknown metabolite.");
                
                unsigned i = found->second;
                NT val = (NT)it.value().get<double>();
                eq_triplets.push_back(Triplet(i, j, val));
            }
        }
        ++j;
    }

    MT A_eq(m, n);
    A_eq.setFromTriplets(eq_triplets.begin(), eq_triplets.end());
    A_eq.makeCompressed();

    VT b_eq = VT::Zero(m);

    // No inequalities besides box bounds
    MT A_in(0, n);
    VT b_in(0);

    return ConstrainedPolytope<Point>(n, A_eq, b_eq, A_in, b_in, lb, ub);
}

// Parses the BiGG JSON model from a given file path.
// @tparam Point the Point type used by the metabolic network
// @param model_path path to the .json file
// @return the parsed metabolic network as a ConstrainedPolytope
template <typename Point>
ConstrainedPolytope<Point> parse_from_json(std::string const& model_path) 
{
    std::ifstream f(model_path);
    
    if (!f.is_open())
        throw std::runtime_error("Cannot open the BiGG model "+model_path);

    nlohmann::json jsn;
    f >> jsn;

    return construct_from_json<Point>(jsn);
}
#endif