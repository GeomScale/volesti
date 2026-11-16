#pragma once

#include "core_types.hpp"
#include <chrono>
#include <string>
#include <vector>

//In this file you can find the definition of supportive functions.
//Timer is a class to count the time for each random walk.
//write to file is used to write the results to a file.
//vector to eigen converts a vector in an eigen compatible form. (eigen is numpy for c++) 
//determine auto walk is used to choose a walk based on the polytope structure.


class Timer {
public:
    Timer(const std::string& name = "");
    void start();
    double stop(const std::string& label = "");
    double get_total_time() const;

private:
    std::string walk_name;
    std::chrono::steady_clock::time_point start_time;
    double total_time;
    bool is_running;
};


void write_to_file(std::string filename, std::vector<Point> const& randPoints);


template <typename MT>
MT vector_to_eigen(const std::vector<Point>& someSamples) {
    if (someSamples.empty()) {
        return MT(); // Return empty matrix if no samples
    }
    
    MT samples(someSamples[0].dimension(), someSamples.size());
    for (size_t jj = 0; jj < someSamples.size(); ++jj) {
        samples.col(jj) = someSamples[jj].getCoefficients();
    }
    return samples;
}

std::string determine_auto_walk(unsigned int dim);