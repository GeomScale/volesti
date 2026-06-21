#pragma once

#include "core_types.hpp"
#include <vector>
#include <string>

// In this file you can find the structures and functions responsible for
// post-processing and summarizing the results of random walk sampling.

struct WalkResult {
    std::vector<Point> samples;
    unsigned int final_ess;
    double generation_time;
    double ess_time;
};


struct WalkStatistics {
    unsigned int final_ess;
    double max_psrf;
    double ks_statistic;
    double ks_p_value;
    double total_time;
    double ess_time;
};

WalkStatistics process_and_print_results(
    const std::vector<Point>& samples, 
    HPOLYTOPE& polytope, 
    const std::string& walk_name, 
    double total_generation_time,
    unsigned int precalculated_ess,
    double total_ess_time
);