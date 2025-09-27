// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and modified by Huu Phuoc Le as part of Google Summer of Code 2022 program
// Contributed and/or modified by Angelos Korakitis, as part of Google Summer of Code 2025 program.


// Licensed under GNU LGPL.3, see LICENCE file


#ifndef VOLESTI_EIGENVALUESPROBLEMS_H
#define VOLESTI_EIGENVALUESPROBLEMS_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

// Spectra library for eigenvalue problems
#include <Spectra/include/Spectra/GenEigsSolver.h>
#include <Spectra/include/Spectra/SymEigsSolver.h>
#include <Spectra/include/Spectra/MatOp/DenseGenMatProd.h>
#include <Spectra/include/Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/include/Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/include/Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/include/Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/include/Spectra/Util/SelectionRule.h>
#include <Spectra/include/Spectra/MatOp/DenseCholesky.h>

#include <Spectra/include/Spectra/MatOp/DenseSymShiftSolve.h>
#include <Spectra/include/Spectra/SymEigsShiftSolver.h>

#include <algorithm>
#include <limits>
#include <cmath>

template<typename NT,
         typename MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>,
         typename VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>>
class EigenvaluesProblems {
private:
    using SparseMT = Eigen::SparseMatrix<NT>;
    
    static constexpr NT eps() { return std::numeric_limits<NT>::epsilon(); }
    static constexpr NT sqrt_eps() { return std::sqrt(eps()); }
    static constexpr NT LARGE_VAL() { return std::numeric_limits<NT>::max(); }
    static constexpr NT EPS() { return std::numeric_limits<NT>::epsilon(); }

    /**
     * @brief Operator for the structured QEP (C1^{-1}C0)
     * Where:
     * C0 = [B2   0]  and C1 = [I   0]
     *      [0    I]           [B1  B0]
     */
    class QEPOperator {
private:
    const MT& A_;  // Coefficient of t^2
    const MT& B_;  // Coefficient of t
    const MT& C_;  // Constant term
    const int n_;
    mutable Eigen::LLT<MT> cholC_;  // For solving C_ systems
    mutable bool chol_computed_ = false;
public:
    using Scalar = NT;

    QEPOperator(const MT& A, const MT& B, const MT& C)
        : A_(A), B_(B), C_(C), n_(A.rows()) {}

    int rows() const { return 2 * n_; }
    int cols() const { return 2 * n_; }

    // Computes y = M^{-1} x, where M = [0 I; -C -B] and x = [x1; A x2]
    void perform_op(const NT* x_in, NT* y_out) const {
    Eigen::Map<const VT> x(x_in, 2 * n_);
    Eigen::Map<VT> y(y_out, 2 * n_);

    // Solve [0, I; -C, -B] * [y1; y2] = [x1; x2]
    // From first row:  0*y1 + I*y2 = x1  =>  y2 = x1
    // From second row: -C*y1 - B*y2 = x2  =>  -C*y1 = x2 + B*y2
    
    y.tail(n_) = x.head(n_);  // y2 = x1
    VT rhs = x.tail(n_) + B_ * y.tail(n_);  // x2 + B*x1

    // Compute Cholesky decomposition if not already done
    if (!chol_computed_) {
        cholC_.compute(C_);
        chol_computed_ = true;
    }

    // Solve C*y1 = -rhs for y1
    if (cholC_.info() == Eigen::Success) {
        y.head(n_) = cholC_.solve(-rhs);
    } else {
        // Fallback to QR decomposition if Cholesky fails
        y.head(n_) = C_.colPivHouseholderQr().solve(-rhs);
    }
}
};

    /**
 * @brief Solve QEP using dense companion matrix approach with improved preconditioning
 * Uses Jacobi preconditioning to improve numerical stability
 */
static NT solveDenseCompanion(const MT& A, const MT& B, const MT& C, VT& eigvec) {
    const int n = A.rows();
    
    try {
        // Build standard companion matrices for QEP: (λ²A + λB + C)v = 0
        MT X(2*n, 2*n), Y(2*n, 2*n);
        X.setZero();
        Y.setZero();
        X.block(0, n, n, n).setIdentity();
        X.block(n, 0, n, n) = -C;
        X.block(n, n, n, n) = -B;
        Y.block(0, 0, n, n).setIdentity();
        Y.block(n, n, n, n) = A;
        
        // Improved preconditioning:
        NT matrix_scale = std::max({A.norm(), B.norm(), C.norm()});
        NT diag_threshold = std::max(NT(1e-14), matrix_scale * eps() * NT(100));
        // Use Jacobi preconditioner
        MT P_inv = Y.diagonal().array().max(diag_threshold).cwiseInverse().matrix().asDiagonal();
        
        // Apply preconditioning: solve P^(-1)*X*u = λ*P^(-1)*Y*u
        // This has the same eigenvalues as the original problem
        X = P_inv * X;
        Y = P_inv * Y;
        
        // Solve generalized eigenvalue problem: X*u = λ*Y*u
        Eigen::GeneralizedEigenSolver<MT> solver(X, Y);
        if (solver.info() != Eigen::Success) {
            return std::numeric_limits<NT>::infinity();
        }
        
        // Find smallest positive real eigenvalue
        NT tolerance = std::max(NT(1e-12), matrix_scale * eps() * NT(100));
        constexpr NT epsBeta = NT(1e-12);
        NT best_lambda = std::numeric_limits<NT>::max();
        int best_idx = -1;
        
        for (int i = 0; i < solver.alphas().size(); ++i) {
            if (std::abs(solver.betas()(i)) < epsBeta) continue;
            
            std::complex<NT> lambda = solver.alphas()(i) / solver.betas()(i);
            if (std::abs(lambda.imag()) > tolerance) continue;
            
            NT real_lambda = lambda.real();
            if (real_lambda > tolerance && real_lambda < best_lambda) {
                best_lambda = real_lambda;
                best_idx = i;
            }
        }
        
        if (best_idx >= 0) {
            // Extract eigenvector (first n components for standard companion form)
            eigvec = solver.eigenvectors().col(best_idx).head(n).real();
            if (eigvec.norm() > eps()) {
                eigvec.normalize();
                return best_lambda;
            }
        }
        
    } catch (const std::exception&) {
        // Handle numerical exceptions gracefully
    }
    
    return std::numeric_limits<NT>::infinity();
}

    /**
     * @brief Solve QEP using shift-invert with dense Spectra solvers only
     */
    static NT solveWithShiftInvert(const MT& A, const MT& B, const MT& C, VT& eigvec) {
        const int n = A.rows();
        
        // Build companion matrices
        MT X(2*n, 2*n), Y(2*n, 2*n);
        X.setZero();
        Y.setZero();
        X.block(0, n, n, n).setIdentity();
        X.block(n, 0, n, n) = -C;
        X.block(n, n, n, n) = -B;
        Y.block(0, 0, n, n).setIdentity();
        Y.block(n, n, n, n) = A;
        
        // Simple shift selection - use a small positive value
        NT matrix_scale = std::max({A.norm(), B.norm(), C.norm()});
        NT shift = std::max(NT(1e-6), matrix_scale * NT(1e-8));
        
        try {
            const int nev = std::min(6, 2*n-2);
            const int ncv = std::min(std::max(2*nev, 20), 2*n-1);
            
            // Use dense shift-invert solver
            Spectra::DenseGenRealShiftSolve<NT> op(X);
            Spectra::GenEigsRealShiftSolver<Spectra::DenseGenRealShiftSolve<NT>> eigs(op, nev, ncv, shift);
            
            eigs.init();
            int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 1000, NT(1e-10));
            
            if (nconv > 0) {
                auto eigenvals = eigs.eigenvalues();
                auto eigenvecs = eigs.eigenvectors();
                
                NT epsPos = std::max(NT(1e-8), matrix_scale * std::numeric_limits<NT>::epsilon() * NT(50));
                NT best = std::numeric_limits<NT>::infinity();
                int best_idx = -1;
                
                for (int i = 0; i < nconv; ++i) {
                    if (std::abs(eigenvals(i).imag()) > epsPos) continue;
                    
                    NT lambda = eigenvals(i).real();
                    if (lambda > epsPos && lambda < best) {
                        best = lambda;
                        best_idx = i;
                    }
                }
                
                if (best_idx >= 0) {
                    eigvec = eigenvecs.col(best_idx).head(n).real();
                    eigvec.normalize();
                    return best;
                }
            }
        } catch (...) {
            return std::numeric_limits<NT>::infinity();
        }
        
        return std::numeric_limits<NT>::infinity();
    }

    /**
     * @brief Solve the structured QEP using Spectra (for larger systems)
     */
    static NT solveQEPSpectra(const MT& B0, const MT& B1, const MT& B2, VT& eigvec) {
        const int m = B0.rows();
        if (m == 0) return NT(0);

        try {
            QEPOperator op(B0, B1, B2);
            const int nev = 6;          // Target 6 eigenvalues
            const int ncv = nev + 3;    // Arnoldi subspace size

            Spectra::GenEigsSolver<QEPOperator> solver(op, nev, ncv);
            solver.init();
            int nconv = solver.compute(Spectra::SortRule::LargestMagn, 2000, NT(1e-10));
            
            if (nconv > 0 && solver.info() == Spectra::CompInfo::Successful) {
                auto eigenvals = solver.eigenvalues();
                auto eigenvecs = solver.eigenvectors();
                
                NT best_t = 0;
                int best_idx = -1;
                const NT tol = NT(1e-8) * B0.norm();
                
                for (int i = 0; i < nconv; ++i) {
                    NT lambda = eigenvals(i).real();  // We want real eigenvalues
                    
                    // We solve for 1/t, so t = 1/lambda
                    if (std::abs(lambda) > tol) {
                        NT t = NT(1)/lambda;
                        if (t > tol && (best_idx == -1 || t < best_t)) {
                            VT candidate = eigenvecs.col(i).real();
                            if (candidate.norm() > sqrt_eps()) {
                                best_t = t;
                                best_idx = i;
                                eigvec = candidate;
                            }
                        }
                    }
                }
                
                if (best_idx >= 0) {
                    eigvec.normalize();
                    return best_t;
                }
            }
        } catch (...) {
            return std::numeric_limits<NT>::infinity();
        }
        
        return std::numeric_limits<NT>::infinity();
    }

    /**
     * @brief Solve the structured QEP with automatic method selection
     */
    static NT solveQEP(const MT& B0, const MT& B1, const MT& B2, VT& eigvec) 
    {
        const int n = B0.rows();
        if (n == 0) return NT(0);
        
        // Use dense companion method for small systems (typically faster and more reliable)
        if (n <= 15) {  // Adjust based on empirical performance (...)
            NT result = solveDenseCompanion(B0, B1, B2, eigvec);
            if (result > NT(0)) {
                return result;
            }
        }
        // else if (n <= 20)
        // {
        //     // For medium-sized systems, try dense Spectra method first
        //     return solveWithShiftInvert(B0, B1, B2, eigvec);
        // }
        return solveQEPSpectra(B0, B1, B2, eigvec);
        
        // Fall back to Spectra-based method for larger systems or if dense method fails
    }

public:
    /**
     * @brief Solve the structured quadratic eigenvalue problem
     * Finds the minimum positive eigenvalue for the problem defined by B0, B1, B2
     * Returns infinity if no finite minimum positive eigenvalue exists or solver fails
     */
    NT minPosQuadraticEigenvalue(const MT& A, const MT& B, const MT& C,
                                 MT& /*X*/, MT& /*Y*/, VT& eigvec,
                                 bool /*updateOnly*/ = false,
                                 int  /*num_constraints*/ = 0)
    {
        int n = A.rows();
        if (n == 0 || A.cols()!=n || B.rows()!=n || B.cols()!=n ||
            C.rows()!=n || C.cols()!=n)
            return std::numeric_limits<NT>::infinity();

        // Test
        // std::cout << "[SOLVER] minPosQuadraticEigenvalue called with size " << n << std::endl;

        NT result = solveQEP(A, B, C, eigvec);
        
        // Ensure we return a valid positive eigenvalue or infinity
        if (result <= NT(0) || !std::isfinite(result)) {
            return std::numeric_limits<NT>::infinity();
        }

        
        return result;
    }

    // /**
    //  * @brief Solve generalized eigenvalue problem B*v = λ*(-A)*v
    //  * Returns minimum positive and maximum negative eigenvalues
    //  */
    // std::pair<NT, NT> symGeneralizedProblem(const MT& A, const MT& B) {
    //     if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
    //         return {std::numeric_limits<NT>::max(), -std::numeric_limits<NT>::max()};
    //     }

    //     // Test
    //     // std::cout << "[SOLVER] symGeneralizedProblem called with size " << A.rows() << std::endl;

    //     try {
    //         // Solve B*v = λ*(-A)*v, which is equivalent to (-A)⁻¹*B*v = (1/λ)*v
    //         MT neg_A = -A;
    //         Eigen::GeneralizedSelfAdjointEigenSolver<MT> solver(B, neg_A);
            
    //         if (solver.info() != Eigen::Success) {
    //             return {std::numeric_limits<NT>::max(), -std::numeric_limits<NT>::max()};
    //         }
            
    //         NT min_pos = std::numeric_limits<NT>::max();
    //         NT max_neg = -std::numeric_limits<NT>::max();
            
    //         auto eigenvals = solver.eigenvalues();
    //         NT tolerance = std::max(NT(1e-12), std::max(A.norm(), B.norm()) * eps() * NT(100000));
            
    //         for (int i = 0; i < eigenvals.size(); ++i) {
    //             NT mu = eigenvals(i);
    //             if (std::abs(mu) < tolerance) continue;
                

    //             NT lambda = 1.0 / mu;
                
    //             if (lambda > tolerance && lambda < min_pos) min_pos = lambda;
    //             if (lambda < -tolerance && lambda > max_neg) max_neg = lambda;
    //         }
            
    //         return {min_pos, max_neg};
            
    //     } catch (...) {
    //         return {std::numeric_limits<NT>::max(), -std::numeric_limits<NT>::max()};
    //     }
    // }


    static std::pair<NT, NT> symGeneralizedProblem(const MT& A, const MT& B) {
        // Check matrix dimensions
        if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
            return {LARGE_VAL(), -LARGE_VAL()};
        }
        
        const int n = A.rows();
        const NT tolerance = std::max(NT(1e-12), 
                 std::max(A.norm(), B.norm()) * eps() * NT(100000));
        
        NT min_pos = LARGE_VAL();
        NT max_neg = -LARGE_VAL();
        
        // First attempt: Use Spectra's iterative solver
        if (n >= 50){
        try {
            const int num_eigs = std::min(10, n);  // Number of eigenvalues to compute
            const int ncv = std::min(2 * num_eigs, n);  // Size of Arnoldi factorization
            
            // For generalized eigenvalue problem Ax = λBx, we solve (A - σB)^(-1)Bx = μx
            // where μ = 1/(λ - σ) and σ is the shift (typically 0)
            MT neg_A = -A;
            
            // Create the shift-solve operator for (A - σB)^(-1) where σ = 0
            // So we need to solve A^(-1)B or in this case (-A)^(-1)B
            Spectra::DenseSymShiftSolve<NT> op(neg_A);  // Only pass the matrix A
            
            // Create the SymEigsShiftSolver with only one template parameter
            Spectra::SymEigsShiftSolver<Spectra::DenseSymShiftSolve<NT>> 
                eigs(op, num_eigs, ncv, 0.0);
            
            eigs.init();
            int nconv = eigs.compute(Spectra::SortRule::LargestMagn);
            
            // Process eigenvalues if Spectra succeeded
            if (eigs.info() == Spectra::CompInfo::Successful && nconv > 0) {
                auto mu_vals = eigs.eigenvalues();
                
                for (int i = 0; i < nconv; ++i) {
                    NT mu = mu_vals(i);
                    // Skip near-zero eigenvalues
                    if (std::abs(mu) < tolerance) continue;
                    
                    // Convert to original eigenvalue: λ = 1/μ
                    NT lambda = NT(1.0) / mu;
                    
                    if (lambda > tolerance && lambda < min_pos) min_pos = lambda;
                    if (lambda < -tolerance && lambda > max_neg) max_neg = lambda;
                }
                
                // Return if we found valid eigenvalues
                if (min_pos != LARGE_VAL() && max_neg != -LARGE_VAL()) {
                    return {min_pos, max_neg};
                }
            }
        } catch (...) {
            // Spectra failed - fall through to Eigen solver
        }
    }
        // Fallback: Use Eigen's dense solver (full eigenvalue decomposition)
        try {
            MT neg_A = -A;
            Eigen::GeneralizedSelfAdjointEigenSolver<MT> solver(B, neg_A);
            
            if (solver.info() != Eigen::Success) {
                return {LARGE_VAL(), -LARGE_VAL()};
            }
            
            auto eigenvals = solver.eigenvalues();
            for (int i = 0; i < eigenvals.size(); ++i) {
                NT lambda = eigenvals(i);
                if (std::abs(lambda) < tolerance) continue;
                if (lambda > tolerance && lambda < min_pos) min_pos = lambda;
                if (lambda < -tolerance && lambda > max_neg) max_neg = lambda;
            }
            
            return {min_pos, max_neg};
        } catch (...) {
            return {LARGE_VAL(), -LARGE_VAL()};
        }
    }



    /**
     * @brief Find minimum eigenvalue of symmetric matrix with robust numerical handling
     * Returns a numerically meaningful minimum eigenvalue for optimization contexts
     */
// NT findSymEigenvalue(const MT& M) {
//     if (M.rows() != M.cols() || M.rows() == 0) {
//         return std::numeric_limits<NT>::infinity();
//     }

//     // Test
//     // std::cout << "[SOLVER] findSymEigenvalue called with size " << M.rows() << std::endl;

//         try {
//             Eigen::SelfAdjointEigenSolver<MT> solver(M, Eigen::EigenvaluesOnly);
//             if (solver.info() != Eigen::Success) {
//                 return std::numeric_limits<NT>::infinity();
//             }
//             return solver.eigenvalues().maxCoeff();
//         } catch (...) {}
        
    
//     // Complete failure
//     return std::numeric_limits<NT>::infinity();
// }


NT findSymEigenvalue(const MT& M) {
    if (M.rows() != M.cols() || M.rows() == 0) {
        return std::numeric_limits<NT>::infinity();
    }
    
    try {
        // Create matrix operation object
        Spectra::DenseSymMatProd<NT> op(M);
        
        // Construct eigen solver object, requesting 1 eigenvalue
        // Parameters: matrix operation, number of eigenvalues, convergence parameter
        // Rule of thumb: convergence parameter = 2 * number of eigenvalues + 1
        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<NT>> solver(op, 1, 7);
        
        // Initialize and compute
        solver.init();
        int nconv = solver.compute(Spectra::SortRule::LargestAlge);
        
        // Check convergence
        if (solver.info() != Spectra::CompInfo::Successful) {
            return std::numeric_limits<NT>::infinity();
        }
        
        // Return the largest eigenvalue
        return solver.eigenvalues()(0);
        
    } catch (...) {
        return std::numeric_limits<NT>::infinity();
    }
}

};

#endif // VOLESTI_EIGENVALUESPROBLEMS_H
