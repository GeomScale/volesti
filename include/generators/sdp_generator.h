// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SDP_GENERATOR_H
#define VOLESTI_SDP_GENERATOR_H

#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iomanip>
// boost
#include <boost/random/normal_distribution.hpp>
// eigen
#include <Eigen/Dense>

typedef boost::mt19937 RNGType;

/// Output format for SDP instance
enum class SDPFormat {
    DENSE,   // Dense format: all elements written
    SPARSE   // SDPA sparse format: only non-zero elements
};

/// Generates a random matrix
template <class NT>
void randomMatrixGOE(Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>& M) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    unsigned m = M.rows();
    boost::normal_distribution<> rdist(0,1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    for (unsigned int i=0; i<m; i++) {
        for (unsigned int j=0; j<m; j++) {
            M(i,j) = rdist(rng);
        }
    }
}

/// Generates a random spectrahedron S(n, m) and saves to file
/// @param filename Output filename
/// @param n Number of decision variables
/// @param m Size of the matrix (m x m)
/// @param format Output format (DENSE or SPARSE)
/// @param tol Tolerance for sparse format (elements with |value| < tol are considered zero)
template<typename NT, typename SpectrahedronType>
void generate_sdp_instance(const std::string& filename, int n, int m, 
                          SDPFormat format = SDPFormat::DENSE, NT tol = 1e-10) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    // Generate random matrices
    MT ones = MT::Ones(m, m);
    MT M = 2 * MT::Random(m, m) - ones;
    MT I = MT::Identity(m, m);
    
    std::vector<MT> matrices(n + 1);
    matrices[0] = -(M * M.transpose()) - I;

    MT ones2 = MT::Ones(m/2, m/2);
    MT MM(m/2, m/2), MMM(m/2, m/2);

    for (int i = 1; i <= n; i++) {
        MM = 2 * MT::Random(m/2, m/2) - ones2;
        MMM.setZero(m/2, m/2);
        
        for (int l = 0; l < m/2; ++l) {
            MMM(l, l) = MM(l, l);
        }
        
        for (int j = 0; j < m/2; ++j) {
            for (int k = j + 1; k < m/2; ++k) {
                MMM(j, k) = MM(j, k);
                MMM(k, j) = MMM(j, k);
            }
        }

        MT A = MT::Zero(m, m);
        for (int j = 0; j < m/2; ++j) {
            for (int k = 0; k < m/2; ++k) {
                A(j, k) = MMM(j, k);
                A(j + m/2, k + m/2) = -MMM(j, k);
            }
        }
        matrices[i] = A;
    }

    // Create LMI and Spectrahedron
    LMI<NT, MT, VT> lmi(matrices);
    SpectrahedronType spectrahedron(lmi);
    
    // Generate random objective function
    VT objective_vec(n);
    for (int i = 0; i < n; i++) {
        objective_vec(i) = (NT)(rand() % 100 - 50) / 10.0;
    }

    // Save to file
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to create SDP instance file: " + filename);
    }
    
    // Set precision for output
    ofs << std::setprecision(6) << std::scientific;
    
    if (format == SDPFormat::SPARSE) {
        // SDPA sparse format
        
        // Line 2: n (number of decision variables)
        ofs << n << std::endl;
        
        // Line 3: number of blocks (always 1)
        ofs << "1" << std::endl;
        
        // Line 4: block size (positive for SDP block)
        ofs << m << std::endl;
        
        // Line 5: objective function coefficients
        ofs << std::fixed;
        for (int i = 0; i < n; i++) {
            ofs << objective_vec(i);
            if (i < n - 1) ofs << " ";
        }
        ofs << std::endl;
        
        // Lines 6+: Sparse matrix entries in format: matno blkno i j value
        // matno: 0 for constant matrix, 1..n for constraint matrices
        // blkno: block number (1-indexed, always 1 in our case)
        // i, j: row and column (1-indexed)
        // Only upper triangular part for symmetric matrices
        ofs << std::scientific;
        
       
for (int mat_idx = 0; mat_idx <= n; mat_idx++) {
    const MT& mat = matrices[mat_idx];
    for (int i = 0; i < m; i++) {
        for (int j = i; j < m; j++) {
            NT value = mat(i, j);
            
            // Sign convention: negate matrices 1..n for SDPA format
            if (mat_idx > 0) {
                value = -value;
            }
            
            if (std::abs(value) >= tol) {
                ofs << mat_idx << " 1 " << (i+1) << " " << (j+1) << " " << value << std::endl;
            }
        }
    }
}
        
        std::cout << "Generated SDP instance (SDPA sparse format): n=" << n << ", m=" << m << std::endl;
        
    } else {
        // Dense format (original)
        ofs << std::fixed;
        
        // Line 1: n (number of decision variables)
        ofs << n << std::endl;
        
        // Line 2: number of blocks (always 1)
        ofs << "1" << std::endl;
        
        // Line 3: m (matrix size)
        ofs << m << std::endl;
        
        // Line 4: objective function coefficients
        for (int i = 0; i < n; i++) {
            ofs << objective_vec(i);
            if (i < n - 1) ofs << " ";
        }
        ofs << std::endl;
        
        // Lines 5+: Write all n+1 matrices in dense format
        for (int mat_idx = 0; mat_idx <= n; mat_idx++) {
            const MT& mat = matrices[mat_idx];
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    // Write each element
                    if (std::abs(mat(i, j)) < tol) {
                        ofs << std::setw(12) << "0";
                    } else {
                        ofs << std::setw(12) << mat(i, j);
                    }
                    if (j < m - 1) ofs << " ";
                }
                ofs << std::endl;
            }
        }
        
        std::cout << "Generated SDP instance (dense format): n=" << n << ", m=" << m << std::endl;
    }
    
    ofs.close();
}

/// Legacy function - generates and returns a Spectrahedron (for in-memory use)
template<typename NT, typename SpectrahedronType>
SpectrahedronType generateSDP(int n, int m) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    MT ones = MT::Ones(m, m);
    MT M = 2 * MT::Random(m, m) - ones;
    MT I = MT::Identity(m, m);
    
    std::vector<MT> matrices(n + 1);
    matrices[0] = -(M * M.transpose()) - I;

    MT ones2 = MT::Ones(m/2, m/2);
    MT MM(m/2, m/2), MMM(m/2, m/2);

    for (int i = 1; i <= n; i++) {
        MM = 2 * MT::Random(m/2, m/2) - ones2;
        MMM.setZero(m/2, m/2);
        
        for (int j = 0; j < m/2; ++j) {
            for (int k = 0; k < m/2; ++k) {
                MMM(j, k) = MM(j, k) + MM(k, j);
            }
        }

        MT A = MT::Zero(m, m);
        for (int j = 0; j < m/2; ++j) {
            for (int k = 0; k < m/2; ++k) {
                A(j, k) = MMM(j, k);
                A(j + m/2, k + m/2) = -MMM(j, k);
            }
        }
        matrices[i] = A;
    }

    LMI<NT, MT, VT> lmi(matrices);
    SpectrahedronType spectrahedron(lmi);
    return spectrahedron;
}

#endif //VOLESTI_SDP_GENERATOR_H