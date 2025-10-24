// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Based on Sieve-SDP facial reduction algorithm

#ifndef VOLESTI_FACIAL_REDUCTION_HPP
#define VOLESTI_FACIAL_REDUCTION_HPP

#include <Eigen/Dense>
#include <vector>
#include <iostream>

/// LMI structure for facial reduction
template<typename NT>
struct FacialReductionLMI {
    using Matrix = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    
    int m;                      ///< Matrix dimension
    int n;                      ///< Number of variables
    Matrix A0;                  ///< Constant matrix
    std::vector<Matrix> A;      ///< Coefficient matrices
    
    FacialReductionLMI() : m(0), n(0) {}
    FacialReductionLMI(int m_, int n_) : m(m_), n(n_), A0(m_, m_), A(n_, Matrix(m_, m_)) {
        A0.setZero();
        for (auto& mat : A) mat.setZero();
    }
};

struct FacialReductionOptions {
    double tolerance = 1e-10;
    int max_iterations = 100;
    bool verbose = false;
};

template<typename NT>
struct FacialReductionResult {
    enum Status { SUCCESS, NO_REDUCTION, INFEASIBLE };
    
    Status status;
    FacialReductionLMI<NT> reduced_lmi;
    int reduced_size;
    std::vector<int> kept_indices;
    std::string message;
    
    FacialReductionResult() : status(NO_REDUCTION), reduced_size(0) {}
};

/// Facial Reduction using Sieve-SDP algorithm
template<typename NT>
class FacialReduction {
public:
    using Matrix = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
    
    explicit FacialReduction(const FacialReductionOptions& opts = FacialReductionOptions())
        : options_(opts) {}
    
    /// Apply facial reduction to LMI
    FacialReductionResult<NT> reduce(const FacialReductionLMI<NT>& lmi) {
        FacialReductionResult<NT> result;
        
        if (options_.verbose) {
            std::cout << "[Facial] Input: m=" << lmi.m << ", n=" << lmi.n << "\n";
        }
        
        std::vector<int> current_indices(lmi.m);
        for (int i = 0; i < lmi.m; ++i) current_indices[i] = i;
        
        FacialReductionLMI<NT> current_lmi = lmi;
        bool reduced = false;
        
        for (int iter = 0; iter < options_.max_iterations; ++iter) {
            auto step_result = reduce_step(current_lmi);
            
            if (step_result.status == FacialReductionResult<NT>::INFEASIBLE) {
                result.status = FacialReductionResult<NT>::INFEASIBLE;
                result.message = "Problem is infeasible";
                return result;
            }
            
            if (step_result.status == FacialReductionResult<NT>::NO_REDUCTION) {
                break;
            }
            
            current_lmi = step_result.reduced_lmi;
            std::vector<int> new_indices;
            for (int idx : step_result.kept_indices) {
                new_indices.push_back(current_indices[idx]);
            }
            current_indices = new_indices;
            reduced = true;
            
            if (options_.verbose) {
                std::cout << "[Facial] Iteration " << iter + 1 
                          << ": reduced to m=" << current_lmi.m << "\n";
            }
        }
        
        if (reduced) {
            result.status = FacialReductionResult<NT>::SUCCESS;
            result.reduced_lmi = current_lmi;
            result.reduced_size = current_lmi.m;
            result.kept_indices = current_indices;
            result.message = "Reduced from " + std::to_string(lmi.m) + " to " + 
                           std::to_string(current_lmi.m) + " dimensions";
        } else {
            result.status = FacialReductionResult<NT>::NO_REDUCTION;
            result.reduced_lmi = lmi;
            result.reduced_size = lmi.m;
            result.kept_indices.resize(lmi.m);
            for (int i = 0; i < lmi.m; ++i) result.kept_indices[i] = i;
            result.message = "No reduction possible";
        }
        
        return result;
    }

private:
    FacialReductionOptions options_;
    
    /// Single reduction step using Sieve-SDP
    FacialReductionResult<NT> reduce_step(const FacialReductionLMI<NT>& lmi) {
        FacialReductionResult<NT> result;
        result.status = FacialReductionResult<NT>::NO_REDUCTION;
        
        // Check A0 and each Ai for block structure
        for (int k = 0; k <= lmi.n; ++k) {
            const Matrix& Ak = (k == 0) ? lmi.A0 : lmi.A[k-1];
            
            // Compute eigendecomposition
            Eigen::SelfAdjointEigenSolver<Matrix> es(Ak);
            const Vector& eigenvals = es.eigenvalues();
            const Matrix& eigenvecs = es.eigenvectors();
            
            // Find negative eigenvalues (< -tolerance)
            std::vector<int> neg_indices;
            for (int i = 0; i < eigenvals.size(); ++i) {
                if (eigenvals[i] < -options_.tolerance) {
                    neg_indices.push_back(i);
                }
            }
            
            // Skip if no structure found
            if (neg_indices.empty() || neg_indices.size() == lmi.m) {
                continue;
            }
            
            // Project to negative eigenspace
            Matrix Q = eigenvecs(Eigen::all, neg_indices);
            int new_m = neg_indices.size();
            
            FacialReductionLMI<NT> reduced(new_m, lmi.n);
            
            reduced.A0 = Q.transpose() * lmi.A0 * Q;
            for (int i = 0; i < lmi.n; ++i) {
                reduced.A[i] = Q.transpose() * lmi.A[i] * Q;
            }
            
            result.status = FacialReductionResult<NT>::SUCCESS;
            result.reduced_lmi = reduced;
            result.reduced_size = new_m;
            result.kept_indices = neg_indices;
            return result;
        }
        
        return result;
    }
};

#endif // VOLESTI_FACIAL_REDUCTION_HPP