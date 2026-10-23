#include "../include/walk_result.hpp"
#include "../include/diagnostics.hpp" 

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <filesystem>

WalkStatistics process_and_print_results(
    const std::vector<Point>& samples, 
    HPOLYTOPE& polytope, 
    const std::string& walk_name, 
    double total_generation_time,
    unsigned int precalculated_ess,
    double total_ess_time,
    unsigned int walk_len,
    const std::string& polytope_name,
    bool show_console_logs) 
{
    // Calculate mixing ratio (Total Steps / ESS)
    double mixing_ratio = 0.0;
    if (precalculated_ess > 0) {
        // total steps = number of saved samples * thinning factor (walk_len) (no burn in included)
        mixing_ratio = static_cast<double>(samples.size() * walk_len) / precalculated_ess;
    }

    // Calculate metrics
    WalkStatistics stats = {precalculated_ess, 0.0, -1.0, -1.0, total_generation_time, total_ess_time, mixing_ratio};

    if (samples.empty()) {
        std::cerr << "[" << walk_name << "] Error: No samples to process!\n";
        return stats;
    }

    //std::cout << "\n--- Processing Statistics for " << walk_name << " ---\n";

    // Convert to Eigen Matrix
    MT samples_mat = vector_to_eigen<MT>(samples);

    std::cout << std::fixed << std::setprecision(4);

    if (show_console_logs) {
    // ESS and mixing rate
        std::cout << "[" << walk_name << "] Final ESS: " << stats.final_ess << "\n";
        std::cout << "[" << walk_name << "] Mixing Ratio (Steps/ESS): " << stats.mixing_ratio << "\n";

        // Time
        std::cout << "[" << walk_name << "] Total Algorithm Time: " << stats.total_time << " seconds\n";
        std::cout << "[" << walk_name << "] Total ESS Time: " << stats.ess_time << " seconds\n";
    }
    // PSRF 

    stats.max_psrf = compute_psrf<NT, VT, MT>(samples);

    if (show_console_logs) {
        std::cout << "[" << walk_name << "] Max PSRF: " << stats.max_psrf << "\n";
    }

    // Condition for KS Test
    // Check if "Gaussian" is in the walk name because KS test is only for uniform
    bool is_gaussian = (walk_name.find("Gaussian") != std::string::npos);

    if (!is_gaussian) {
        auto ks_results = compute_ks_test<HPOLYTOPE, MT>(polytope, samples_mat, stats.final_ess);
        stats.ks_statistic = ks_results.ks_stat;
        stats.ks_p_value = ks_results.p_val;

        if (show_console_logs) {
            std::cout << "[" << walk_name << "] KS Statistic: " << stats.ks_statistic << "\n";
            std::cout << "[" << walk_name << "] P-Value:      " << stats.ks_p_value << "\n";
        }
    } else {
        if (show_console_logs) {
            std::cout << "[" << walk_name << "] KS Test:      Skipped (Gaussian Distribution)\n";
        }
    }
    std::cout << "--------------------------------------------------\n";

    // Append results to the benchmark CSV file 
    std::string filename = "benchmark_results.csv";
    
    // Check if file exists before we open it in append mode
    bool file_exists = std::filesystem::exists(filename);

    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app); // Append mode
    
    if (outfile.is_open()) {
        // If this is the very first time creating the file, write the header row
        if (!file_exists) {
            outfile << "Polytope,Dimension,Method,Time_Sec,Points,ESS,Mixing_Ratio,Max_PSRF,KS_Stat,KS_P_Value\n";
        }

        // Write the data row
        outfile << polytope_name << ", "
                << polytope.dimension() << ", " 
                << walk_name << ", "
                << std::fixed << std::setprecision(4) << total_generation_time << ", "
                << samples.size() << ", "
                << stats.final_ess << ", "
                << std::fixed << std::setprecision(4) << stats.mixing_ratio << ", " 
                << std::fixed << std::setprecision(4) << stats.max_psrf << ", "
                << std::fixed << std::setprecision(4) << stats.ks_statistic << ", "
                << std::fixed << std::setprecision(4) << stats.ks_p_value << "\n";
                
        outfile.close();
    } else {
        std::cerr << "!!! Unable to open " << filename << " to save data.\n";
    }

    return stats;
}