// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

// The test cases are taken from https://github.com/unc-optimization/SieveSDP

#include "doctest.h"
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <cmath>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"
#include "optimization/simulated_annealing.hpp"
#include "preprocess/spectrahedron/facial_reduction.hpp"

using NT = double;
using Kernel = Cartesian<NT>;
using Point = point<Kernel>;
using SPECTRAHEDRON = Spectrahedron<Point>;

enum class FRTestResult {
    PASS,
    FAIL,
    REDUCTION_FAILED
};

struct FRTestCase {
    std::string filename;
    NT expected_minimum;
    NT relative_error;
    NT tolerance;
    std::string description;
};

struct SDPSolution {
    NT objective_value;
    double solve_time_ms;
    int matrix_size;
    int num_variables;
};

FacialReductionLMI<NT> extract_lmi_from_spectrahedron(const SPECTRAHEDRON& spec) {
    const auto& volesti_lmi = spec.getLMI();
    int d = volesti_lmi.dimension();
    int m = volesti_lmi.sizeOfMatrices();
    const auto& all_matrices = volesti_lmi.getMatrices();
    
    FacialReductionLMI<NT> fr_lmi(m, d);
    fr_lmi.A0 = all_matrices.empty() ? Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>::Zero(m, m) : all_matrices[0];
    
    for (int i = 0; i < d; ++i) {
        fr_lmi.A[i] = (i + 1 < all_matrices.size()) ? all_matrices[i + 1] : 
                      Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>::Zero(m, m);
    }
    
    return fr_lmi;
}

SPECTRAHEDRON build_spectrahedron_from_lmi(const FacialReductionLMI<NT>& fr_lmi) {
    using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
    using LMI_TYPE = LMI<NT, MT, VT>;
    
    std::vector<MT> all_matrices;
    all_matrices.reserve(fr_lmi.n + 1);
    all_matrices.push_back(fr_lmi.A0);
    for (int i = 0; i < fr_lmi.n; ++i) {
        all_matrices.push_back(fr_lmi.A[i]);
    }
    
    return SPECTRAHEDRON(LMI_TYPE(all_matrices));
}

SPECTRAHEDRON load_problem(const std::string& filename, Point& objective) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    SPECTRAHEDRON spec;
    SdpaFormatManager<NT> sdpa_manager;
    sdpa_manager.loadSDPAFormatFile(file, spec, objective);
    file.close();
    
    return spec;
}

SDPSolution solve_problem(SPECTRAHEDRON& spec, const Point& objective, NT relative_error) {
    auto start = std::chrono::high_resolution_clock::now();
    
    Point initial_point(spec.getLMI().dimension());
    SimulatedAnnealingSettings<Point> settings(relative_error);
    Point solution;
    
    NT result = solve_sdp(spec, objective, settings, initial_point, solution, false);
    
    auto end = std::chrono::high_resolution_clock::now();
    double time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    return {result, time_ms, static_cast<int>(spec.getLMI().sizeOfMatrices()), static_cast<int>(spec.getLMI().dimension())};
}

bool apply_facial_reduction(const SPECTRAHEDRON& input_spec, const Point& input_obj,
                           SPECTRAHEDRON& output_spec, Point& output_obj,
                           int& original_size, int& reduced_size) {
    try {
        FacialReductionLMI<NT> fr_lmi = extract_lmi_from_spectrahedron(input_spec);
        original_size = fr_lmi.m;
        
        FacialReductionOptions options;
        options.verbose = false;
        options.tolerance = 1e-8;
        
        FacialReduction<NT> reducer(options);
        FacialReductionResult<NT> result = reducer.reduce(fr_lmi);
        
        if (result.status != FacialReductionResult<NT>::Status::SUCCESS) {
            return false;
        }
        
        reduced_size = result.reduced_size;
        output_spec = build_spectrahedron_from_lmi(result.reduced_lmi);
        output_obj = input_obj;
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

static FRTestResult run_facial_reduction_test(const FRTestCase& test_case) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "\n=== Testing " << test_case.description << " ===" << std::endl;
    std::cout << "File: " << test_case.filename << std::endl;
    
    try {
        Point objective_original, objective_reduced;
        
        // Solve original problem
        std::cout << "\n[Original Problem]" << std::endl;
        SPECTRAHEDRON spec_original = load_problem(test_case.filename, objective_original);
        SDPSolution original = solve_problem(spec_original, objective_original, test_case.relative_error);
        
        std::cout << "  Matrix size: " << original.matrix_size << " x " << original.matrix_size << std::endl;
        std::cout << "  Variables: " << original.num_variables << std::endl;
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  Objective: " << original.objective_value << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Time: " << original.solve_time_ms / 1000.0 << " s" << std::endl;
        
        // Apply facial reduction
        std::cout << "\n[Facial Reduction]" << std::endl;
        SPECTRAHEDRON spec_input = load_problem(test_case.filename, objective_reduced);
        
        auto start_reduction = std::chrono::high_resolution_clock::now();
        SPECTRAHEDRON spec_reduced;
        Point obj_reduced;
        int original_size, reduced_size;
        
        if (!apply_facial_reduction(spec_input, objective_reduced, spec_reduced, obj_reduced, 
                                   original_size, reduced_size)) {
            std::cerr << "  ERROR: Facial reduction failed" << std::endl;
            return FRTestResult::REDUCTION_FAILED;
        }
        
        auto end_reduction = std::chrono::high_resolution_clock::now();
        double reduction_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_reduction - start_reduction).count();
        
        std::cout << "  Original: " << original_size << " x " << original_size << std::endl;
        std::cout << "  Reduced: " << reduced_size << " x " << reduced_size << std::endl;
        std::cout << "  Eliminated: " << (original_size - reduced_size) << " (" 
                  << std::fixed << std::setprecision(1)
                  << (100.0 * (original_size - reduced_size) / original_size) << "%)" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Time: " << reduction_time_ms / 1000.0 << " s" << std::endl;
        
        // Solve reduced problem
        std::cout << "\n[Reduced Problem]" << std::endl;
        SDPSolution reduced = solve_problem(spec_reduced, obj_reduced, test_case.relative_error);
        
        std::cout << "  Matrix size: " << reduced.matrix_size << " x " << reduced.matrix_size << std::endl;
        std::cout << "  Variables: " << reduced.num_variables << std::endl;
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  Objective: " << reduced.objective_value << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Time: " << reduced.solve_time_ms / 1000.0 << " s" << std::endl;
        
        // Compare results
        const NT difference = std::abs(original.objective_value - reduced.objective_value);
        const double total_reduced_time = (reduction_time_ms + reduced.solve_time_ms) / 1000.0;
        const double speedup = original.solve_time_ms / (reduction_time_ms + reduced.solve_time_ms);
        
        std::cout << "\n[Comparison]" << std::endl;
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "  Expected: " << test_case.expected_minimum << std::endl;
        std::cout << "  Original: " << original.objective_value 
                  << " (error: " << std::abs(original.objective_value - test_case.expected_minimum) << ")" << std::endl;
        std::cout << "  Reduced:  " << reduced.objective_value 
                  << " (error: " << std::abs(reduced.objective_value - test_case.expected_minimum) << ")" << std::endl;
        std::cout << "  Difference: " << difference << std::endl;
        std::cout << "  Tolerance: " << test_case.tolerance << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Total time (reduced): " << total_reduced_time << " s" << std::endl;
        std::cout << "  Speedup: " << speedup << "x" << std::endl;
        
        // Determine result
        FRTestResult result;
        std::string result_string;
        
        if (difference <= test_case.tolerance) {
            result = FRTestResult::PASS;
            result_string = "PASS";
        } else {
            result = FRTestResult::FAIL;
            result_string = "FAIL";
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "\nTest result: " << result_string << std::endl;
        std::cout << "Total test time: " << total_duration.count() / 1000.0 << " seconds" << std::endl;
        std::cout << std::string(60, '-') << std::endl;
        
        return result;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return FRTestResult::FAIL;
    }
}

static void test_fr_3x5() {
    FRTestCase test_case = {
        .filename = "./dataset/example_1.txt",
        .expected_minimum = 0,
        .relative_error = 1e-2,
        .tolerance =1e-1,
        .description = "3x5 SDP with facial reduction"
    };
    
    FRTestResult result = run_facial_reduction_test(test_case);
    CHECK(result == FRTestResult::PASS);
    if (result == FRTestResult::REDUCTION_FAILED) {
        WARN("Facial reduction failed for this instance");
    }
}

static void test_fr_2x3() {
    FRTestCase test_case = {
        .filename = "./dataset/example_2.txt",
        .expected_minimum = -1.9710734588e+00,
        .relative_error = 1e-3,
        .tolerance = 1e-1,
        .description = "2x3 SDP with facial reduction"
    };
    
    FRTestResult result = run_facial_reduction_test(test_case);
    CHECK(result == FRTestResult::PASS);
    if (result == FRTestResult::REDUCTION_FAILED) {
        WARN("Facial reduction failed for this instance");
    }
}

static void test_fr_3x3() {
    FRTestCase test_case = {
        .filename = "./dataset/example_4.txt",
        .expected_minimum = -0,
        .relative_error = 1e-3,
        .tolerance = 1e-1,
        .description = "3x3 SDP with facial reduction"
    };
    
    FRTestResult result = run_facial_reduction_test(test_case);
    CHECK(result == FRTestResult::PASS);
    if (result == FRTestResult::REDUCTION_FAILED) {
        WARN("Facial reduction failed for this instance");
    }
}

static void test_fr_20x20() {
    FRTestCase test_case = {
        .filename = "./dataset/example_3.txt",
        .expected_minimum = -0,
        .relative_error = 1e-3,
        .tolerance = 1e-1,
        .description = "20x20 SDP with facial reduction"
    };
    
    FRTestResult result = run_facial_reduction_test(test_case);
    CHECK(result == FRTestResult::PASS);
    if (result == FRTestResult::REDUCTION_FAILED) {
        WARN("Facial reduction failed for this instance");
    }
}


TEST_CASE("Facial Reduction Preprocessing") {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Testing Facial Reduction Preprocessing" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    SUBCASE("(3x5)") {
        test_fr_3x5();
    }

    SUBCASE("(2x3)") {
        test_fr_2x3();
    }

    SUBCASE("(3x3)") {
        test_fr_3x3();
    }

    SUBCASE("(20x20)") {
        test_fr_20x20();
    }
}