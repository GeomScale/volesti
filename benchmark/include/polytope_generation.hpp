#pragma once

#include "core_types.hpp"
#include "walk_parameters.hpp"

#include <string>
#include <map>
#include <functional>
#include <iostream>
#include <filesystem>
#include <stdexcept>

#include "known_polytope_generators.h"
#include "order_polytope_generator.h"
#include "custom_generators.h"

using PolytopeGeneratorFn = std::function<HPOLYTOPE(unsigned int, const BenchmarkConfig&)>;

inline HPOLYTOPE create_polytope(const std::string& choice, unsigned int dim, const BenchmarkConfig& config) {

    static std::map<std::string, PolytopeGeneratorFn> factory = {
        {"Cube", [](unsigned int d, const BenchmarkConfig&) {
            return generate_cube<HPOLYTOPE>(d, false);
        }},
        {"Simplex", [](unsigned int d, const BenchmarkConfig&) {
            return generate_simplex<HPOLYTOPE>(d, false);
        }},
        {"Birkhoff", [](unsigned int d, const BenchmarkConfig&) {
            return generate_birkhoff<HPOLYTOPE>(d);
        }},
        {"Cross", [](unsigned int d, const BenchmarkConfig&) {
            return generate_cross<HPOLYTOPE>(d, false);
        }},
        {"OrderPolytope", [](unsigned int d, const BenchmarkConfig& cfg) {
            unsigned int m = 3 * d;
            int seed = cfg.base_seed + d;
            return random_orderpoly<HPOLYTOPE, double>(d, m, seed);
        }},
{       "Custom", [](unsigned int /*ignored_d*/, const BenchmarkConfig& cfg) {
            if (cfg.custom_A_file.empty() || cfg.custom_b_file.empty()) {
                throw std::runtime_error("\n[CRITICAL ERROR] Custom polytope chosen but CSV file paths are missing in the config JSON!\n");
            }

            // verify files exist on disk before touching them
            if (!std::filesystem::exists(cfg.custom_A_file)) {
                throw std::runtime_error(
                    "\n[CRITICAL ERROR] Cannot find A matrix file: " + cfg.custom_A_file + 
                    "\n-> Hint: Ensure the CSV is in the same directory you are running the executable from, or use an absolute path in your config.\n"
                );
            }
            if (!std::filesystem::exists(cfg.custom_b_file)) {
                throw std::runtime_error(
                    "\n[CRITICAL ERROR] Cannot find b vector file: " + cfg.custom_b_file + 
                    "\n-> Hint: Ensure the CSV is in the same directory you are running the executable from, or use an absolute path in your config.\n"
                );
            }

            return load_custom_polytope<HPOLYTOPE>(cfg.custom_A_file, cfg.custom_b_file);
        }}
    };

    auto it = factory.find(choice);
    if (it != factory.end()) {
        return it->second(dim, config); // Call the matched lambda
    } else {
        std::cerr << "!!! Unknown polytope choice: " << choice << ". Defaulting to Cube.\n";
        return generate_cube<HPOLYTOPE>(dim, false);
    }
}