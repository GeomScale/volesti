#pragma once

#include <string>
#include <map>
#include <vector>

// In this file you can find functions for configuration structures and helper
// functions used to control benchmark experiments for random walks.

// A struct to hold the specific settings for a single random walk
struct WalkSettings {
    bool enabled = true;
    unsigned int samples = 1000;
    unsigned int walk_len_multiplier = 0;
    unsigned int walk_len_base = 1;

    // For Gaussian
    double a_i_param = 1.0;
};

// The global configuration object
struct BenchmarkConfig {
    
    unsigned int target_ESS = 500;
    double time_limit_sec = 1200.0;
    int base_seed = 42;
    unsigned int dimension = 3;
    std::vector<unsigned int> dimensions;
    double angle = 0.0;
    std::string polytope_choice = "Cube";
    std::string custom_A_file;
    std::string custom_b_file;
    bool use_dynamic_batch = true;
    bool write_to_file = false;
    bool rounding = false;
    std::string rounding_method = "max_ellipsoid";
    bool auto_walk = false;
    bool show_console_logs = false;
    bool show_menu = true;

    std::map<std::string, WalkSettings> walk_settings;
};



// Function declaration to load the JSON file
BenchmarkConfig load_benchmark_config(const std::string& filepath);

// Helper functions that the runners will use
unsigned int get_initial_batch_size(const std::string& walk_name, const BenchmarkConfig& config);
unsigned int compute_dynamic_walk_len(const std::string& walk_name, unsigned int dim, const BenchmarkConfig& config);