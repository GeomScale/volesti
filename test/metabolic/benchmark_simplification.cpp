// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "io/bigg_parser.hpp"
#include "preprocess/metabolic/simplification_exhaustive.hpp"
#include "preprocess/metabolic/simplification_clarkson.hpp"
#include "lp_oracles/metabolic_polyoracles.hpp"
#include <algorithm>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <vector>

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;

unsigned bounds_relaxed(Polytope const& Po, Polytope const& Ps) {
    return Po.getNumFiniteBounds()-Ps.getNumFiniteBounds();
}

unsigned dims_fixed(Polytope const& Po, Polytope const& Ps) {
    return Ps.getNumEqualities()-Po.getNumEqualities();
}

void print_row(std::string method,
               Polytope const& Po,
               Polytope const& Ps,
               std::string status,
               double elapsed_s)
{
    std::cout << " " << std::left << std::setw(12) << method
              << std::right << std::setw(8) << bounds_relaxed(Po, Ps)
              << std::setw(8) << dims_fixed(Po, Ps)
              << std::setw(10) << status
              << std::setw(10) << std::fixed << std::setprecision(3)
              << elapsed_s << "s" << std::endl;
}

template <typename Simplifier, typename Config>
void run_method(std::string method, Polytope const& P, bool dimension_fixing) {
    Config config;
    config.fix_dimensions = dimension_fixing;

    auto start = std::chrono::high_resolution_clock::now();
    Simplifier simplifier(P, config);
    auto [Ps, success] = simplifier.simplify();
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_s = std::chrono::duration<double>(end-start).count();

    auto equal_orac = are_equal(P, Ps);

    std::string status = equal_orac.value ? "OK"
                       : !equal_orac.solved ? "UNVALIDATED"
                                            : "NOT EQUAL";

    print_row(method, P, Ps, status, elapsed_s);
}

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

    run_method<ExhaustiveSimplifier<Point>, ExhaustiveConfig>(
        "exhaustive", P, dimension_fixing);
    
    run_method<ClarksonSimplifier<Point>, ClarksonConfig>(
        "clarkson", P, dimension_fixing);
}

int main() {
    // Collects and sorts the models to maintain a fixed order in testing.
    std::vector<std::filesystem::path> models;
    for (auto const& file : std::filesystem::directory_iterator(BIGG_DIR)) {
        if (file.path().extension() != ".json") continue;
        models.push_back(file.path());
    }
    std::sort(models.begin(), models.end());

    for (auto const& model : models) {
        // benchmark(model.string(), false); // without dimension fixing
        benchmark(model.string(), true);  // with dimension fixing
    }

    return 0;
}