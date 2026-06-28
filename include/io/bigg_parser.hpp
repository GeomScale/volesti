// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BIGG_PARSER_HPP
#define BIGG_PARSER_HPP

#include "convex_bodies/metabolic_polytope.h"
#include "io/json.hpp"
#include <unordered_map>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <string>
#include <limits>

namespace bigg {
    // Parses a BiGG JSON model into a MetabolicPolytope.
    // The model is represented as:
    // reactions (n) - variables,
    // metabolites (m) - equality rows,
    // reaction lower/upper bound - b_l/b_u vectors.
    // @tparam Point the Point type used by MetabolicPolytope
    // @param jsn a parsed nlohmann::json object holding the model
    // @return the MetabolicPolytope
    template <typename Point>
    MetabolicPolytope<Point> construct_from_json(nlohmann::json const& jsn) {
        typedef MetabolicPolytope<Point> Polytope;
        typedef typename Polytope::MT MT;
        typedef typename Polytope::VT VT;
        typedef typename Polytope::NT NT;
        typedef typename Polytope::Triplet Triplet;

        const NT INF = std::numeric_limits<NT>::infinity();

        if (!jsn.contains("metabolites") || !jsn.contains("reactions"))
            throw std::runtime_error("Not BiGG file: missing `metabolites` or `reactions` field");
            
        auto const& metabolites = jsn.at("metabolites");
        auto const& reactions = jsn.at("reactions");
        unsigned m = (unsigned)metabolites.size(); // metabolite count
        unsigned n = (unsigned)reactions.size();   // reaction count

        MT A_eq(m, n);
        VT b_eq = VT::Zero(m);
        VT b_l(n), b_u(n);

        // Maps the metabolites to the integers in [0,m)
        std::unordered_map<std::string, unsigned> metabolite_index;
        {
            unsigned i = 0;
            for (auto const& metabolite : metabolites) {
                std::string met_id = metabolite.at("id").get<std::string>();
                metabolite_index.emplace(met_id, i++);
            }
        }

        // Builds the stoichiometric matrix A_eq and the vectors b_l, b_u
        std::vector<Triplet> triplets;
        unsigned j = 0;
        for (auto const& reaction : reactions) {
            NT low = reaction.contains("lower_bound") ? reaction.at("lower_bound").get<NT>() : -INF;
            NT high = reaction.contains("upper_bound") ? reaction.at("upper_bound").get<NT>() : INF;
        
            b_l(j) = low;
            b_u(j) = high;

            if (reaction.contains("metabolites")) {
                for (auto it = reaction.at("metabolites").begin(); it != reaction.at("metabolites").end(); ++it) {
                    auto found = metabolite_index.find(it.key());
                    if (found == metabolite_index.end()) // metabolite wasn't found, file is corrupt
                        throw std::runtime_error("Reaction references unknown metabolite.");
                    
                    unsigned i = found->second;
                    NT val = (NT)it.value().get<double>();
                    triplets.push_back(Triplet(i, j, val));
                }
            }
            ++j;
        }

        A_eq.setFromTriplets(triplets.begin(), triplets.end());
        A_eq.makeCompressed();

        return Polytope(n, A_eq, b_l, b_u, b_eq);
    }

    // Parses the BiGG JSON model from a given file path.
    // @tparam Point the Point type used by MetabolicPolytope
    // @param model_path path to the .json file
    // @return the parsed MetabolicPolytope
    template <typename Point>
    MetabolicPolytope<Point> parse_from_json(std::string const& model_path) {
        std::ifstream f(model_path);
        
        if (!f.is_open())
            throw std::runtime_error("Cannot open the BiGG model "+model_path);

        nlohmann::json jsn;
        f >> jsn;

        return construct_from_json<Point>(jsn);
    }
}
#endif