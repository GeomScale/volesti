#pragma once

#include <iostream>
#include <string>
#include <iomanip> // Required for std::setprecision

inline void draw_progress_bar(const std::string& walk_name, 
                              unsigned int current, 
                              unsigned int total, 
                              unsigned int total_samples_so_far,
                              unsigned int current_ESS,
                              unsigned int walk_len,
                              int bar_width = 25) 
{
    if (total == 0) return; // Prevent division by zero

    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(bar_width * progress);
    
    // Clear the line and start drawing
    std::cout << "\r" << std::string(100, ' ') << "\r[" << walk_name << "] Batch: [";
    
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    
    // Print percentages, fractions
    std::cout << "] " << static_cast<int>(progress * 100.0) << "% (" 
              << current << "/" << total << ")";

    // Calculate and print the live mixing ratio
    if (current_ESS > 0) {
        double live_mixing_ratio = static_cast<double>(total_samples_so_far * walk_len) / current_ESS;
        std::cout << " | Mix Ratio: " << std::fixed << std::setprecision(2) << live_mixing_ratio;
    } else {
        std::cout << " | Mix Ratio: N/A"; // Display N/A during the very first batch before ESS is known
    }

    std::cout << std::flush;
}