#pragma once

#include <map>
#include <string>
#include <functional>

#include "core_types.hpp"
#include "walk_run.hpp"
#include "walk_result.hpp"

/*
 * This file allows sampling methods to be selected by name at runtime,
 * instead of creating specific walk types directly in the code.
 *
 * Every walk is registered with a string name and exposed through the same
 * function interface. This makes it possible to store and execute different
 * walk implementations in a uniform way.
 *
 * RunFunction defines the common function signature that all registered walks
 * must follow.
 *
 * WalkRegistry stores the mapping between walk names and their corresponding
 * execution functions. It acts as a lookup table, making it easy to choose a
 * sampling method based on user input or benchmark settings.
 *
 * The execute_walk() helper adapts template-based walk implementations so they
 * can be stored in the registry and called through the common interface.
 *
 * Main functions:
 *
 *   - register_walk(): Registers a walk in the registry.
 *   - get_walk_registry(): Returns the global registry.
 *   - initialize_all_walks(): Registers all available walks.
 *
 * Together with walk_adapters.hpp, this file provides the infrastructure for
 * selecting and running sampling methods through a single, consistent interface.
 */

// Each walk must match this callable signature
using RunFunction = std::function<WalkResult(
    HPOLYTOPE&,
    const Point&,
    RNGType&,
    const BenchmarkConfig&,
    const std::string&
)>;

// Registry type alias
using WalkRegistry = std::map<std::string, RunFunction>;

// Function that returns the global registry
WalkRegistry& get_walk_registry();

// Helper to register a walk. It binds the name to the specific function.
void register_walk(const std::string& name, RunFunction fn);

// Executes the walk by calling the sample_using_walk
template <typename WalkType>
WalkResult execute_walk(
    HPOLYTOPE& P,
    const Point& c,
    RNGType& r,
    const BenchmarkConfig& cfg,
    const std::string& name
) {
    return sample_using_walk<WalkType>(P, c, r, cfg, name);
}

void initialize_all_walks();