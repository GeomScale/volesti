// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"
#include "optimization/simulated_annealing.hpp"

using NT = double;
using Kernel = Cartesian<NT>;
using Point = point<Kernel>;
using SPECTRAHEDRON = Spectrahedron<Point>;

/**
 * Test result enumeration
 */
enum class SDPTestResult {
    PASS,
    FAIL,
    INFEASIBLE
};

/**
 * Test case configuration for SDP solver tests
 */
struct SDPTestCase {
    std::string filename;
    NT expected_minimum;
    NT relative_error;
    NT additive_error_factor;
    std::string description;
};

/**
 * Generic function to test SDP solver on a given test case
 * @param test_case Configuration for the test
 * @return SDPTestResult indicating the test outcome
 */
static SDPTestResult run_sdp_test(const SDPTestCase& test_case) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "\n=== Testing " << test_case.description << " ===" << std::endl;
    std::cout << "File: " << test_case.filename << std::endl;

    // Load SDP problem instance
    std::ifstream input_file(test_case.filename);
    if (!input_file.is_open()) {
        std::cerr << "ERROR: Failed to open SDP instance file: " << test_case.filename << std::endl;
        REQUIRE_MESSAGE(false, "Failed to open SDP instance file");
        return SDPTestResult::FAIL;
    }

    SPECTRAHEDRON spectrahedron;
    Point objective_function;
    SdpaFormatManager<NT> sdpa_manager;
    
    try {
        sdpa_manager.loadSDPAFormatFile(input_file, spectrahedron, objective_function);
    } catch (const std::exception& e) {
        std::cerr << "Error loading SDPA file: " << e.what() << std::endl;
        return SDPTestResult::FAIL;
    }

    // Initialize solver parameters
    Point initial_point(spectrahedron.getLMI().dimension());  // Origin as starting point
    SimulatedAnnealingSettings<Point> solver_settings(test_case.relative_error);

    // Calculate tolerance for success check
    const NT absolute_tolerance = test_case.additive_error_factor * std::abs(test_case.expected_minimum);
    
    // Use scientific notation for output
    std::cout << std::scientific << std::setprecision(16);
    std::cout << "Expected minimum: " << test_case.expected_minimum << std::endl;

    // Solve the SDP problem
    Point solution;
    NT computed_minimum = solve_sdp(spectrahedron,
                                   objective_function,
                                   solver_settings,
                                   initial_point,
                                   solution,
                                   /*verbose=*/ false); // change to true for detailed output

    // Evaluate results
    const NT absolute_error = std::abs(test_case.expected_minimum - computed_minimum);
    
    SDPTestResult result;
    std::string result_string;
    
    if (computed_minimum < test_case.expected_minimum) {
        result = SDPTestResult::INFEASIBLE;
        result_string = "INFEASIBLE";
    } else if (absolute_error <= absolute_tolerance) {
        result = SDPTestResult::PASS;
        result_string = "PASS";
    } else {
        result = SDPTestResult::FAIL;
        result_string = "FAIL";
    }

    std::cout << "Computed minimum: " << computed_minimum << std::endl;
    std::cout << "Tolerance: ±" << absolute_tolerance << std::endl;
    std::cout << "Absolute error: " << absolute_error << std::endl;
    std::cout << "Test result: " << result_string << std::endl;
    
    // Calculate and display elapsed time
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Elapsed time: " << duration.count() / 1000.0 << " seconds" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    return result;
}

// Individual test case wrappers for better test organization
static void test_sdp_2x8() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp__2_8_d.txt",
        .expected_minimum = -1.3888e+00,
        .relative_error = 1e-2,
        .additive_error_factor = 5e-2, 
        .description = "2x8 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_20x20() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp_prob_20_20_d.txt",
        .expected_minimum = -1.6535356e+00,
        .relative_error = 1e-3,
        .additive_error_factor = 5e-2,  
        .description = "20x20 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_20x30() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp__20_30_d.txt",
        .expected_minimum = -1.89076e+01,
        .relative_error = 1e-3,
        .additive_error_factor = 5e-2, 
        .description = "20x30 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_50x40() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp__50_40_d.txt",
        .expected_minimum = -3.1513966779338073e+01,
        .relative_error = 1e-3,
        .additive_error_factor = 5e-2, 
        .description = "50x40 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_20x100() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp__20_100_d.txt",
        .expected_minimum = -6.74974e+00,
        .relative_error = 1e-3,
        .additive_error_factor = 5e-2,
        .description = "20x100 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_50x200() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp__50_200_d.txt",
        .expected_minimum = -7.71954e+00,
        .relative_error = 1e-3,
        .additive_error_factor = 5e-2,
        .description = "50x200 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_200x15() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp_prob_200_15_d.txt",
        .expected_minimum = -3.1610166002684946e+12,
        .relative_error = 1e-4,
        .additive_error_factor = 5e-2, 
        .description = "200x15 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_400x20() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp_prob_400_20_d.txt",
        .expected_minimum = -8.7728166762094482e+11,
        .relative_error = 1e-4,
        .additive_error_factor = 5e-2, 
        .description = "400x20 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

static void test_sdp_600x25() {
    SDPTestCase test_case = {
        .filename = "../data/spectra_data/sdp_prob_600_25_d.txt",
        .expected_minimum = -1.3547280395306393e+11,
        .relative_error = 1e-2,
        .additive_error_factor = 5e-2,  
        .description = "600x25 SDP instance"
    };
    
    SDPTestResult result = run_sdp_test(test_case);
    CHECK(result != SDPTestResult::FAIL);
    if (result == SDPTestResult::INFEASIBLE) {
        WARN("Test detected infeasible solution (computed minimum < expected minimum)");
    }
}

/**
 * Main test suite for SDP solver using simulated annealing
 * Tests various problem sizes to validate solver robustness
 */
TEST_CASE("SDP Solver - Simulated Annealing") {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Testing SDP Solver with Simulated Annealing" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    SUBCASE("(2x8)") {
        test_sdp_2x8();
    }

    SUBCASE("(20x20)") {
        test_sdp_20x20();
    }

    SUBCASE("(20x30)") {
        test_sdp_20x30();
    }

    SUBCASE("(50x40)") {
        test_sdp_50x40();
    }

    SUBCASE("(20x100)") {
        test_sdp_20x100();
    }

    SUBCASE("(50x200)") {
        test_sdp_50x200();
    }

    SUBCASE("(200x15)") {
        test_sdp_200x15();
    }

    SUBCASE("(400x20)") {
        test_sdp_400x20();
    }

    SUBCASE("(600x25)") {
        test_sdp_600x25();
    }
}