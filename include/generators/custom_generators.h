#ifndef CUSTOM_GENERATORS_HPP
#define CUSTOM_GENERATORS_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <Eigen/Dense>

// -------------------------------------------------------------------------
// Helper: Reads a CSV file into an Eigen Matrix
// -------------------------------------------------------------------------
template <typename NT>
Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> read_csv_to_eigen(const std::string &path) {
    std::ifstream indata;
    indata.open(path);
    
    if (!indata.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }

    std::string line;
    std::vector<NT> values;
    unsigned int rows = 0;
    
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            // Check for empty cells usually caused by trailing commas
            if (!cell.empty()) {
                values.push_back(static_cast<NT>(std::stod(cell)));
            }
        }
        ++rows;
    }
    
    if (rows == 0) return Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>();

    // Calculate columns
    unsigned int cols = values.size() / rows;
    
    // Map the std::vector to an Eigen Matrix
    // We use RowMajor because CSVs are read row by row
    return Eigen::Map<const Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(values.data(), rows, cols);
}

// -------------------------------------------------------------------------
// Generator: Loads a Polytope defined by Ax <= b from CSV files
// -------------------------------------------------------------------------
template <typename Polytope>
Polytope load_custom_polytope(const std::string &file_A, const std::string &file_b) {
    
    // Define types based on the Polytope template
    typedef typename Polytope::NT NT; // Number Type (e.g., double)
    typedef typename Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT; // Matrix Type
    typedef typename Polytope::VT VT; // Vector Type
    
    std::cout << "Loading polytope from CSVs..." << std::endl;

    // Load matrices using the helper
    MT A_raw = read_csv_to_eigen<NT>(file_A);
    MT b_raw = read_csv_to_eigen<NT>(file_b);
    
    // Safety check
    if (A_raw.rows() != b_raw.rows()) {
        throw std::runtime_error("Dimension mismatch: Rows in A do not match rows in b.");
    }

    // Convert b from Matrix (Nx1) to Vector (N)
    VT b = b_raw.col(0); 
    
    unsigned int dim = A_raw.cols();
    unsigned int num_constraints = A_raw.rows();
    
    std::cout << "Successfully loaded: " << num_constraints << " constraints in " << dim << " dimensions." << std::endl;

    // Return the Polytope
    return Polytope(dim, A_raw, b);
}

#endif // CUSTOM_GENERATORS_HPP