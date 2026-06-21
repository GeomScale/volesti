#pragma once

#include <map>
#include <string>
#include <functional>

#include "core_types.hpp"
#include "walk_run.hpp"
#include "walk_result.hpp"

/*
 * The goal here is to decouple walk selection from walk implementation. Rather than
 * directly instantiating specific walk types throughout the codebase, every
 * sampling method is registered under a string identifier and exposed through
 * a common callable interface.
 *
 * The RunFunction type defines the standardized signature that every walk
 * execution routine must follow. This allows all methods, regardless of their
 * underlying implementation, to be stored and invoked uniformly.
 *
 * The WalkRegistry maps method names to executable functions. The registry 
 * acts as a lookup table that enables runtime selection of sampling algorithms 
 * based on user configuration, command-line arguments, or benchmark settings.
 *
 * The execute_walk() helper wraps the generic sample_using_walk() function,
 * converting a walk type into a registry-compatible callable. This allows
 * template-based walk implementations to be registered and invoked through
 * the same runtime interface.
 *
 * The registration functions provide:
 *
 *   - register_walk(): Adds a walk implementation to the registry.
 *   - get_walk_registry(): Accesses the global registry instance.
 *   - initialize_all_walks(): Registers all available walk methods at startup.
 *
 * Together with walk_adapters.hpp, this file forms the core abstraction layer
 * of the framework: walk_adapters.hpp unifies algorithm interfaces, while
 * walk_registry.hpp enables dynamic selection and execution of any registered
 * sampling method through a common lookup mechanism.
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