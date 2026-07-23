#include "../include/benchmark_utils.hpp"
#include <iostream>
#include <fstream>

PushBackWalkPolicy push_back_policy;

// Timer Class. This is usefull to avoid using new time variables every time we need to count the time.
Timer::Timer(const std::string& name) : walk_name(name), total_time(0.0) {}

void Timer::start() { 
    start_time = std::chrono::steady_clock::now(); 
    is_running = true;
}

double Timer::stop(const std::string& label) {
    auto end_time = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count();
    total_time += elapsed;
    is_running = false;
    
    // Only print if a label was provided
    if (!label.empty()) {
        std::cout << "[" << walk_name << "] " << label << " = " << elapsed << " s\n";
    }
    return elapsed;
}

double Timer::get_total_time() const { 

    if (is_running) {
        auto current_time = std::chrono::steady_clock::now();
        double current_elapsed = std::chrono::duration<double>(current_time - start_time).count();
        return total_time + current_elapsed;
    }
    
    return total_time; 
}

void write_to_file(std::string filename, std::vector<Point> const& randPoints) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }

    // Save current cout buffer and redirect to the file
    auto coutbuf = std::cout.rdbuf(out.rdbuf()); 
    
    for(size_t i = 0; i < randPoints.size(); ++i) {
        randPoints[i].print();
    }
    
    // Reset cout back to standard output
    std::cout.rdbuf(coutbuf); 
}

std::string determine_auto_walk(unsigned int dim) {
    if (dim >= 1 && dim <= 10) {
        return "BallWalk";
    } else if (dim >= 11 && dim <= 20) {
        return "RDHRWalk";
    } else if (dim >= 21 && dim <= 30) {
        return "BilliardWalk";
    } else if (dim >= 31 && dim <= 49) {
        return "CDHRWalk";
    } else { 
        return "AcceleratedBilliardWalk";
    }
}