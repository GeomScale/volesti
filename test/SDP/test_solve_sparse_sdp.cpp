// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "convex_bodies/spectrahedra/sparse_spectrahedron.h"
#include "SDPAFormatManager.h"
#include "optimization/simulated_annealing.hpp"
#include "generators/sdp_generator.h"

using NT = double;
using Kernel = Cartesian<NT>;
using Point = point<Kernel>;
using DENSE_SPECTRAHEDRON = Spectrahedron<Point>;
using SPARSE_SPECTRAHEDRON = SparseSpectrahedron<Point>;

/**
 * Test result with timing information
 */
struct SDPTestResult {
    enum Status { PASS, FAIL, INFEASIBLE };
    
    Status status;
    NT computed_minimum;
    NT absolute_error;
    double elapsed_seconds;
    
    bool passed() const { return status == PASS; }
};

/**
 * Test case configuration
 */
struct SDPTestCase {
    std::string filename;
    NT expected_minimum;
    NT relative_error;
    NT additive_error_factor;
    std::string description;
};

/**
 * Generic SDP test runner - works with BOTH dense and sparse
 */
template <typename SpectrahedronType>
static SDPTestResult run_sdp_test_generic(
    const SDPTestCase& test_case,
    const std::string& implementation_name
) {

    // create an sdp test input and save to file
    // filename format: sdp__m_d.txt
    // if the file already exists, it will not be overwritten
    // Check if file exists, generate if needed
    if (!std::ifstream(test_case.filename).good()) {
        std::cout << "Generating SDP instance file: " << test_case.filename << std::endl;
        
        // Parse m and d from filename (format: sdp__m_d.txt or sdp_prob_m_d.txt)
        size_t pos1 = test_case.filename.find("__");
        if (pos1 == std::string::npos) {
            pos1 = test_case.filename.rfind("_", test_case.filename.rfind("_") - 1);
        }
        size_t pos2 = test_case.filename.find("_", pos1 + 2);
        size_t pos3 = test_case.filename.find(".txt");
        
        if (pos1 == std::string::npos || pos2 == std::string::npos || pos3 == std::string::npos) {
            throw std::runtime_error("Invalid filename format. Expected: sdp__m_d.txt or sdp_prob_m_d.txt");
        }
        
        int m = std::stoi(test_case.filename.substr(pos1 + 2, pos2 - (pos1 + 2)));
        int d = std::stoi(test_case.filename.substr(pos2 + 1, pos3 - (pos2 + 1)));
        
        // Generate and save the SDP instance
        generate_sdp_instance<double, SpectrahedronType>(
        test_case.filename, m, d, SDPFormat::SPARSE);
        std::cout << "SDP instance file generated successfully." << std::endl;
    }


    std::cout << "\n=== " << implementation_name << ": " 
              << test_case.description << " ===" << std::endl;
    std::cout << "File: " << test_case.filename << std::endl;

    // Load problem
    std::ifstream input_file(test_case.filename);
    REQUIRE_MESSAGE(input_file.is_open(), "Failed to open SDP instance file");

    SpectrahedronType spectrahedron;
    Point objective_function;
    SdpaFormatManager<NT> sdpa_manager;
    
    try {
        sdpa_manager.loadSDPAFormatFile(input_file, spectrahedron, objective_function);
    } catch (const std::exception& e) {
        std::cerr << "Error loading SDPA file: " << e.what() << std::endl;
        return {SDPTestResult::FAIL, NT(0), NT(0), 0.0};
    }
    

// Origin as initial 
Point initial_point = Point(spectrahedron.getLMI().dimension());

spectrahedron.set_interior_point(initial_point);

SimulatedAnnealingSettings<Point> solver_settings(test_case.relative_error);
const NT absolute_tolerance = test_case.additive_error_factor * 
                              std::abs(test_case.expected_minimum);

// Solve with timing
auto start = std::chrono::high_resolution_clock::now();

Point solution;
NT computed_minimum;


    try {
        computed_minimum = solve_sdp(
            spectrahedron,
            objective_function,
            solver_settings,
            initial_point,
            solution,
            false  // verbose=false for cleaner test output
        );
    } catch (const std::exception& e) {
        std::cerr << "Error in solve_sdp: " << e.what() << std::endl;
        return {SDPTestResult::FAIL, NT(0), NT(0), 0.0};
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    // Evaluate results
    const NT absolute_error = std::abs(test_case.expected_minimum - computed_minimum);
    
    SDPTestResult::Status status;
    if (computed_minimum < test_case.expected_minimum) {
        status = SDPTestResult::INFEASIBLE;
    } else if (absolute_error <= absolute_tolerance) {
        status = SDPTestResult::PASS;
    } else {
        status = SDPTestResult::FAIL;
    }

    std::cout << "Expected:  " << test_case.expected_minimum << std::endl;
    std::cout << "Computed:  " << computed_minimum << std::endl;
    std::cout << "Error:     " << absolute_error << std::endl;
    std::cout << "Time:      " << elapsed << " seconds" << std::endl;
    std::cout << "Status:    " << (status == SDPTestResult::PASS ? "PASS" : 
                                   status == SDPTestResult::INFEASIBLE ? "INFEASIBLE" : "FAIL") 
              << std::endl;

    return {status, computed_minimum, absolute_error, elapsed};
}

/**
 * Test suite for sparse SDP implementation
 */
TEST_CASE("Sparse SDP Implementation") {

    SUBCASE("(40x300)") {
        SDPTestCase test_case = {
            .filename = "../spectra_data/sdp__40_300_s.txt",
            .expected_minimum = -4.36959,
            .relative_error = 1e-3,
            .additive_error_factor = 5e-2,
            .description = "40x300 Sparse SDP instance"
        };
        auto result = run_sdp_test_generic<SPARSE_SPECTRAHEDRON>(
            test_case, "SPARSE"
        );
        
        CHECK(result.passed());

    }

    SUBCASE("(40x100)") {
        SDPTestCase test_case = {
            .filename = "../spectra_data/sdp__40_100_s.txt",
            .expected_minimum = -14.262,
            .relative_error = 1e-3,
            .additive_error_factor = 5e-2,
            .description = "40x100 Sparse  SDP instance"
        };
        auto result = run_sdp_test_generic<SPARSE_SPECTRAHEDRON>(
            test_case, "SPARSE"
        );
        
        CHECK(result.passed());
    }

    SUBCASE("(30x100)") {
        SDPTestCase test_case = {
            .filename = "../spectra_data/sdp__30_100_s.txt",
            .expected_minimum = -10.4938,
            .relative_error = 1e-3,
            .additive_error_factor = 5e-2,
            .description = "30x100 Sparse SDP instance"
        };
        
        auto result = run_sdp_test_generic<SPARSE_SPECTRAHEDRON>(
            test_case, "SPARSE"
        );
        
        CHECK(result.passed());

    }
    

    SUBCASE("(50x200)") {
        SDPTestCase test_case = {
            .filename = "../spectra_data/sdp__50_200_s.txt",
            .expected_minimum = -9.58509,
            .relative_error = 1e-3,
            .additive_error_factor = 5e-2,
            .description = "50x200 Sparse SDP instance"
        };
        
        auto result = run_sdp_test_generic<SPARSE_SPECTRAHEDRON>(
            test_case, "SPARSE"
        );
        
        CHECK(result.passed());

    }

}   