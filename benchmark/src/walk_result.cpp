#include "../include/walk_result.hpp"
#include "../include/diagnostics.hpp" 

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

WalkStatistics process_and_print_results(
    const std::vector<Point>& samples, 
    HPOLYTOPE& polytope, 
    const std::string& walk_name, 
    double total_generation_time,
    unsigned int precalculated_ess,
    double total_ess_time) 
{
    // Calculate metrics
    WalkStatistics stats = {precalculated_ess, 0.0, 0.0, 0.0, total_generation_time, total_ess_time};

    if (samples.empty()) {
        std::cerr << "[" << walk_name << "] Error: No samples to process!\n";
        return stats;
    }

    //std::cout << "\n--- Processing Statistics for " << walk_name << " ---\n";

    // Convert to Eigen Matrix
    MT samples_mat = vector_to_eigen<MT>(samples);

    // ESS
    std::cout << "[" << walk_name << "] Final ESS: " << stats.final_ess << "\n";

    // Time
    std::cout << "[" << walk_name << "] Total Algorithm Time: " << stats.total_time << " seconds\n";
    std::cout << "[" << walk_name << "] Total ESS Time: " << stats.ess_time << " seconds\n";

    // PSRF 
    stats.max_psrf = compute_psrf<NT, VT, MT>(samples);
    std::cout << "[" << walk_name << "] Max PSRF: " << stats.max_psrf << "\n";

    // KS Test
    auto ks_results = compute_ks_test<HPOLYTOPE, MT>(polytope, samples_mat, stats.final_ess);

    stats.ks_statistic = ks_results.ks_stat;
    stats.ks_p_value = ks_results.p_val;

    std::cout << "[" << walk_name << "] KS Statistic: " << stats.ks_statistic << "\n";
    std::cout << "[" << walk_name << "] P-Value:      " << stats.ks_p_value << "\n";
    std::cout << "--------------------------------------------------\n";

    // Append results to the benchmark CSV file
    std::ofstream outfile;
    outfile.open("benchmark_results.txt", std::ios_base::app); // Append mode
    
    if (outfile.is_open()) {
        // Format: Dim, Method, Time(s), Points, ESS, PSRF, KS_Stat, P_Val
        outfile << polytope.dimension() << ", " 
                << walk_name << ", "
                << std::fixed << std::setprecision(4) << total_generation_time << ", "
                << samples.size() << ", "
                << stats.final_ess << ", "
                << stats.max_psrf << ", "
                << stats.ks_statistic << ", "
                << stats.ks_p_value << "\n";
        outfile.close();
    } else {
        std::cerr << "!!! Unable to open benchmark_results.txt to save data.\n";
    }

    return stats;
}