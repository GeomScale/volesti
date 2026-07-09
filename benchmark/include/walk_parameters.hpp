#pragma once

#include <string>
#include <map>
#include <vector>

// In this file you can find functions for configuration structures and helper
// functions used to control benchmark experiments for random walks.

// A struct to hold the specific settings for a single random walk
struct WalkSettings {
    bool enabled;
    unsigned int samples;
    unsigned int walk_len_multiplier;
    unsigned int walk_len_base;

    // For Gaussian
    double a_i_param = 1.0;
};

// The global configuration object
struct BenchmarkConfig {
    
    unsigned int target_ESS;
    double time_limit_sec;
    int base_seed;
    unsigned int dimension;
    std::vector<unsigned int> dimensions;
    double angle;
    std::string polytope_choice;
    std::string custom_A_file;
    std::string custom_b_file;
    bool use_dynamic_batch;
    bool write_to_file;
    bool rounding;

    std::map<std::string, WalkSettings> walk_settings;
};



// Function declaration to load the JSON file
BenchmarkConfig load_benchmark_config(const std::string& filepath);

// Helper functions that the runners will use
unsigned int get_initial_batch_size(const std::string& walk_name, const BenchmarkConfig& config);
unsigned int compute_dynamic_walk_len(const std::string& walk_name, unsigned int dim, const BenchmarkConfig& config);