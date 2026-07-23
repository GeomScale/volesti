#pragma once

#include <algorithm> 
#include <tuple>
#include <vector>

#include "core_types.hpp" 
#include "benchmark_utils.hpp" 

#include "sampling/sample_correlation_matrices.hpp"
#include "matrix_operations/EigenvaluesProblems.h"
#include "diagnostics/effective_sample_size.hpp"
#include "diagnostics/univariate_psrf.hpp"
#include "diagnostics/scaling_ratio.hpp"
#include "diagnostics/KS_test.hpp"

// In this file you can find utility functions for computing diagnostic
// metrics on samples produced by random walks.

// The file works with Eigen matrices and uses helper utilities such as
// vector_to_eigen to convert sample containers into matrix form.

// Computes ESS on a given Eigen matrix and returns the min ESS
template <typename NT, typename VT, typename MT>
unsigned int compute_ess(const MT& samples) {
    unsigned int min_ess = 0;
    // call internal ESS function
    VT ess_vector = effective_sample_size<NT, VT, MT>(samples, min_ess);
    return min_ess;
}

// Computes PSRF for all accumulated points
template <typename NT, typename VT, typename MT>
double compute_psrf(const std::vector<Point>& someSamples) {
    // Convert the vector of points to an Eigen matrix using our utility
    MT finalSamples = vector_to_eigen<MT>(someSamples);
    
    // Call Volesti's internal PSRF function
    VT psrf = univariate_psrf<NT, VT, MT>(finalSamples);
    return psrf.maxCoeff();
}


// Struct to neatly pass the KS results back
struct KSTestResult {
    double ks_stat;
    double p_val;
    std::vector<double> observed;
    std::vector<double> expected;
};

// Computes the KS Test and scaling ratios
template <typename PolytopeType, typename MT>
KSTestResult compute_ks_test(const PolytopeType& polytope, const MT& samples_mat, double current_ESS) {
    // Calculate a safe thinning factor based on ESS
    double safe_ess = (current_ESS > 0) ? current_ESS : 1.0;
    int computed_thin = static_cast<int>(samples_mat.cols() / safe_ess);
    int thin_factor = std::max(10, computed_thin * 2);
    
    // Call KS test
    auto [ks_stat, p_val, observed, expected] = global_scaling_test(polytope, samples_mat, thin_factor);  

    return {ks_stat, p_val, observed, expected};
}