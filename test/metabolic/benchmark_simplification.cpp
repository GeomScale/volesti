// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "io/bigg_parser.hpp"
#include "preprocess/metabolic/exhaustive_simplification.hpp"
#include "preprocess/metabolic/clarkson_simplification.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <vector>

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;
typedef typename Polytope::MT MT;
typedef typename Polytope::VT VT;

// Prints a row of the benchmark table.
// @param method the name of the method
// @param bounds_relaxed the number of bounds relaxed
// @param dims_fixed the number of dimensions fixed
// @param success whether the simplification was successful
// @param elapsed_ms the time taken in milliseconds
void print_row(char const* method,
               unsigned bounds_relaxed,
               unsigned dims_fixed,
               bool success,
               double elapsed_ms)
{
    std::cout << " " << std::left << std::setw(12) << method
              << std::right << std::setw(8) << bounds_relaxed
              << std::setw(8) << dims_fixed
              << std::setw(10) << (success ? "OK" : "FAIL")
              << std::setw(10) << std::fixed << std::setprecision(3)
              << elapsed_ms << "s" << std::endl;
}

// Runs both simplification methods on the same model and reports the results.
// @param model the path to the JSON file
// @param dimension_fixing whether degenerate dimensions should be fixed
void benchmark(std::string const& model, bool dimension_fixing) {
    Polytope P = parse_from_json<Point>(model);

    std::cout << "\n --- " << std::filesystem::path(model).stem().string()
              << " (n = " << P.getDimension()
              << ", m = " << P.getNumEqualities()
              << ", finite bounds = " << P.getNumFiniteBounds()
              << ", dimension fixing = " << (dimension_fixing ? "true" : "false")
              << ") ---" << std::endl;


    std::cout << " " << std::left << std::setw(12) << "method"
              << std::right << std::setw(8) << "bounds"
              << std::setw(8) << "dims"
              << std::setw(10) << "status"
              << std::setw(11) << "time" << std::endl;

    
    exhaustive_simplification::Config exhaustive_config;
    exhaustive_config.fix_dimensions = dimension_fixing;
    
    auto ex_start = std::chrono::high_resolution_clock::now();
    auto ex_result = exhaustive_simplification::simplify(P, exhaustive_config);
    auto ex_end = std::chrono::high_resolution_clock::now();
    double ex_elapsed_s = std::chrono::duration<double>(ex_end-ex_start).count();

    print_row("exhaustive", ex_result.bounds_relaxed, ex_result.dims_fixed, ex_result.success, ex_elapsed_s);
        
    clarkson_simplification::Config clarkson_config;
    clarkson_config.fix_dimensions = dimension_fixing;

    auto cl_start = std::chrono::high_resolution_clock::now();
    auto cl_result = clarkson_simplification::simplify(P, clarkson_config);
    auto cl_end = std::chrono::high_resolution_clock::now();
    double cl_elapsed_s = std::chrono::duration<double>(cl_end-cl_start).count();
    
    print_row("clarkson", cl_result.bounds_relaxed, cl_result.dims_fixed, cl_result.success, cl_elapsed_s);
}

int main() {
    for (auto const& file : std::filesystem::directory_iterator(BIGG_DIR)) {
        if (file.path().extension() != ".json") continue;
        benchmark(file.path().string(), false); // without dimension fixing
        benchmark(file.path().string(), true);  // with dimension fixing
    }

    return 0;
}