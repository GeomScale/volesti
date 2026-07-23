#pragma once

#include <iostream>
#include <string>
#include <iomanip>

inline void draw_progress_bar(const std::string& walk_name, 
                              unsigned int current, 
                              unsigned int total, 
                              unsigned int total_samples_so_far,
                              unsigned int current_ESS,
                              unsigned int walk_len,
                              double elapsed_seconds,
                              double estimated_remaining_seconds,
                              int bar_width = 25) 
{
    if (total == 0) return; 

    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(bar_width * progress);
    
    // clear the line and start drawing
    std::cout << "\r" << std::string(120, ' ') << "\r[" << walk_name << "] Batch: [";
    //std::cout << "\33[2K\r[" << walk_name << "] Batch: [";
    
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    
    // print percentages, fractions
    std::cout << "] " << static_cast<int>(progress * 100.0) << "% (" 
              << current << "/" << total << ")";

    // calculate and print the live mixing ratio
    if (current_ESS > 0) {
        double live_mixing_ratio = static_cast<double>(total_samples_so_far * walk_len) / current_ESS;
        std::cout << " | Mix Ratio: " << std::fixed << std::setprecision(2) << live_mixing_ratio;
    } else {
        std::cout << " | Mix Ratio: N/A"; // Display N/A during the very first batch before ESS is known
    }

    // Print Elapsed Time 
    unsigned int e_total_secs = static_cast<unsigned int>(elapsed_seconds);
    unsigned int e_hours = e_total_secs / 3600;
    unsigned int e_minutes = (e_total_secs % 3600) / 60;
    unsigned int e_seconds = e_total_secs % 60;

    std::cout << " | Elapsed: ";
    if (e_hours > 0) std::cout << e_hours << "h " << e_minutes << "m";
    else if (e_minutes > 0) std::cout << e_minutes << "m " << e_seconds << "s";
    else std::cout << e_seconds << "s";

    // Print ETA
    if (estimated_remaining_seconds >= 0.0) {
        unsigned int total_secs = static_cast<unsigned int>(estimated_remaining_seconds);
        unsigned int hours = total_secs / 3600;
        unsigned int minutes = (total_secs % 3600) / 60;
        unsigned int seconds = total_secs % 60;

        std::cout << " | ETA: ";
        if (hours > 0) {
            std::cout << hours << "h " << minutes << "m";
        } else if (minutes > 0) {
            std::cout << minutes << "m " << seconds << "s";
        } else {
            std::cout << seconds << "s";
        }
    } else {
        std::cout << " | ETA: Calculating...";
    }

    std::cout << std::flush;
}