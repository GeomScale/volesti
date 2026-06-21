#pragma once

#include <iostream>
#include <string>

inline void draw_progress_bar(const std::string& walk_name, 
                              unsigned int current, 
                              unsigned int total, 
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
    
    // Print percentages, fractions, and flush to console
    std::cout << "] " << static_cast<int>(progress * 100.0) << "% (" 
              << current << "/" << total << ")" << std::flush;
}