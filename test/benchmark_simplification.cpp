// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/metabolic_polytope.h"
#include "io/bigg_parser.hpp"
#include "simplification/warm_start.hpp"
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

void benchmark(std::string const& model, bool dimension_fixing) {
    Polytope P = parse_from_json<Point>(model);

    simplification::Config config;
    config.fix_dimensions = dimension_fixing;

    auto start = std::chrono::high_resolution_clock::now();
    auto result = simplification::simplify(P, config);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end-start).count();

    std::cout << "\n--- "<< model
              << " dimension_fixing = " << (dimension_fixing ? "true" : "false")
              << " ---\n";
    std::cout << " n (reactions)   : " << P.getDimension() << "\n";
    std::cout << " m (metabolites) : " << P.getEqualities().rows() << "\n";
    std::cout << " bounds_relaxed  : " << result.bounds_relaxed << "\n";
    std::cout << " dims_fixed      : " << result.dims_fixed << "\n";
    std::cout << " time            : " << std::fixed 
                                       << std::setprecision(3)
                                       << elapsed << "s\n";
}

int main() {
    for (auto const& model : std::filesystem::directory_iterator(BIGG_DIR)) {
        benchmark(model.path().string(), false); // without dimension fixing
        benchmark(model.path().string(), true);  // with dimension fixing
    }

    return 0;
}