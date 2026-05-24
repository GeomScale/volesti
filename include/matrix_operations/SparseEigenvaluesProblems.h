// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SPARSE_EIGENVALUES_PROBLEMS_H
#define VOLESTI_SPARSE_EIGENVALUES_PROBLEMS_H

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

// Spectra library for sparse eigenvalue problems
#include <Spectra/include/Spectra/GenEigsSolver.h>
#include <Spectra/include/Spectra/SymGEigsSolver.h>
#include <Spectra/include/Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/include/Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/include/Spectra/MatOp/SparseCholesky.h>

/// Solver for various eigenvalue problems arising in convex optimization (Sparse Version)
/// Provides methods for quadratic eigenvalue problems (QEP), generalized eigenvalue
/// problems, and negative definiteness checks via eigenvalue analysis.
/// Uses sparse matrix representations for improved efficiency on large-scale problems.
/// \tparam NT Numeric type for scalar values
/// \tparam SparseMT Sparse matrix type (defaults to Eigen::SparseMatrix)
/// \tparam VT Vector type (defaults to Eigen::Matrix column vector with dynamic size)
template<typename NT,
         typename SparseMT = Eigen::SparseMatrix<NT>,
         typename VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>>
class SparseEigenvaluesProblems {
private:
    /// Machine epsilon for numeric type NT
    static constexpr NT eps() { return std::numeric_limits<NT>::epsilon(); }
    
    /// Square root of machine epsilon (used for moderate tolerance checks)
    static constexpr NT sqrt_eps() { return std::sqrt(eps()); }
    
    /// Largest representable value for numeric type NT
    static constexpr NT LARGE_VAL() { return std::numeric_limits<NT>::max(); }
    
    /// Machine epsilon alias for consistency with existing code
    static constexpr NT EPS() { return std::numeric_limits<NT>::epsilon(); }

    /// Operator for structured quadratic eigenvalue problem (QEP) linearization (Sparse Version)
    /// Implements the matrix-vector product for the linearized system C^{-1}(Bx1 + x2)
    /// where the original QEP is: (A + λB + λ²C)v = 0
    class SparseQEPOperator {
    private:
        const SparseMT& A_;  ///< Coefficient matrix A in QEP (sparse)
        const SparseMT& B_;  ///< Coefficient matrix B in QEP (sparse)
        const SparseMT& C_;  ///< Coefficient matrix C in QEP (must be positive definite, sparse)
        const int n_;  ///< Dimension of original problem
        mutable Eigen::SimplicialLLT<SparseMT> lltC_;    ///< Sparse Cholesky decomposition of C
        mutable Eigen::SimplicialLDLT<SparseMT> ldltC_;  ///< Sparse LDLT decomposition of C (fallback)
        mutable bool chol_computed_ = false;       ///< Flag indicating if decomposition is computed
        mutable bool use_llt_ = false;             ///< Flag indicating which decomposition to use
        
    public:
        using Scalar = NT;  ///< Scalar type required by Spectra library

        /// Construct QEP operator from coefficient matrices
        /// \param[in] A Constant term matrix (sparse)
        /// \param[in] B Linear term matrix (sparse)
        /// \param[in] C Quadratic term matrix (should be positive definite, sparse)
        SparseQEPOperator(const SparseMT& A, const SparseMT& B, const SparseMT& C)
            : A_(A), B_(B), C_(C), n_(A.rows()) {}

        /// Number of rows in linearized system (twice original dimension)
        int rows() const { return 2 * n_; }
        
        /// Number of columns in linearized system (twice original dimension)
        int cols() const { return 2 * n_; }

        /// Perform matrix-vector product y = Op * x for Spectra eigenvalue solver
        /// Implements the linearized QEP system operator using efficient sparse factorizations
        /// \param[in] x_in Input vector of size 2*n
        /// \param[out] y_out Output vector of size 2*n
        void perform_op(const NT* x_in, NT* y_out) const {
            Eigen::Map<const VT> x(x_in, 2*n_);
            Eigen::Map<VT> y(y_out, 2*n_);
            
            // Compute factorization once on first call (lazy initialization)
            // Try Simplicial LLT first (faster for positive definite sparse matrices)
            if (!chol_computed_) {
                lltC_.compute(C_);
                if (lltC_.info() == Eigen::Success) {
                    use_llt_ = true;
                } else {
                    ldltC_.compute(C_);
                    use_llt_ = false;
                }
                chol_computed_ = true;
            }
            
            // Compute right-hand side: rhs = B * x.head(n) + x.tail(n)
            VT rhs(n_);
            rhs = B_ * x.head(n_);
            rhs += x.tail(n_);
            
            // Solve C * y.head(n) = -rhs using precomputed factorization
            if (use_llt_) {
                y.head(n_) = lltC_.solve(-rhs);
            } else if (ldltC_.info() == Eigen::Success) {
                y.head(n_) = ldltC_.solve(-rhs);
            } else {
                // Fallback to direct sparse solve (extremely rare, for robustness)
                Eigen::SparseLU<SparseMT> lu;
                lu.compute(C_);
                if (lu.info() == Eigen::Success) {
                    y.head(n_) = lu.solve(-rhs);
                } else {
                    y.head(n_).setZero();
                }
            }
            
            // Set bottom half of output: y.tail(n) = x.head(n)
            y.tail(n_) = x.head(n_);
        }
    };

    /// Solve quadratic eigenvalue problem to find smallest positive parameter (Sparse Version)
    /// Solves (B0 + t*B1 + t²*B2)v = 0 for the smallest positive t
    /// Uses Spectra library with linearization approach via SparseQEPOperator
    /// \param[in] B0 Constant coefficient matrix (sparse)
    /// \param[in] B1 Linear coefficient matrix (sparse)
    /// \param[in] B2 Quadratic coefficient matrix (sparse)
    /// \param[out] eigvec Eigenvector corresponding to smallest positive t
    /// \return Smallest positive value of t, or infinity if none exists
    static NT solveSparseQEP(const SparseMT& B0, const SparseMT& B1, const SparseMT& B2, VT& eigvec) {
        const int m = B0.rows();
        if (m == 0) return std::numeric_limits<NT>::infinity();
        
        try {
            SparseQEPOperator op(B0, B1, B2);
            
            // Empirical choice: request 2 eigenvalues to improve chances of finding valid positive t
            // for faster convergence use nev = 1, but more unstable, may miss valid eigenvalue
            const int nev = 1; 
            const int ncv = std::min(10, 2*m);
            
            Spectra::GenEigsSolver<SparseQEPOperator> solver(op, nev, ncv);
            solver.init();
            
            int nconv = solver.compute(Spectra::SortRule::LargestMagn, 1000, NT(1e-10));
            
            if (nconv > 0 && solver.info() == Spectra::CompInfo::Successful) {
                auto eigenvals = solver.eigenvalues();
                auto eigenvecs = solver.eigenvectors();
                
                const NT tol = std::max(NT(1e-8) * B0.norm(), eps());
                
                // Find largest positive lambda (smallest positive t)
                NT best_lambda = NT(0);
                int best_idx = -1;
                
                for (int i = 0; i < nconv; ++i) {
                    // Skip complex eigenvalues
                    if (std::abs(eigenvals(i).imag()) > tol) continue;
                    
                    NT lambda = eigenvals(i).real();
                    
                    // We want largest positive lambda (gives smallest t = 1/lambda)
                    if (lambda > tol && lambda > best_lambda) {
                        VT candidate = eigenvecs.col(i).real().head(m);
                        
                        if (candidate.squaredNorm() > tol * tol) {
                            best_lambda = lambda;
                            best_idx = i;
                            eigvec = candidate;
                        }
                    }
                }
                
                if (best_idx >= 0) {
                    eigvec.normalize();
                    return NT(1) / best_lambda;
                }
            }
        } catch (...) {
            return std::numeric_limits<NT>::infinity();
        }
        
        return std::numeric_limits<NT>::infinity();
    }

public:
    /// Find minimum positive eigenvalue for quadratic eigenvalue problem (Sparse Version)
    /// Solves (A + t*B + t²*C)v = 0 to find smallest positive t
    /// Used in convex optimization for barrier function step size computation
    /// \param[in] A Constant coefficient matrix (sparse)
    /// \param[in] B Linear coefficient matrix (sparse)
    /// \param[in] C Quadratic coefficient matrix (sparse)
    /// \param[in,out] X Auxiliary matrix (unused, kept for API compatibility)
    /// \param[in,out] Y Auxiliary matrix (unused, kept for API compatibility)
    /// \param[out] eigvec Eigenvector corresponding to minimum positive eigenvalue
    /// \param[in] updateOnly Unused flag (kept for API compatibility)
    /// \param[in] num_constraints Unused parameter (kept for API compatibility)
    /// \return Minimum positive t satisfying the QEP, or infinity if none exists
    NT minPosQuadraticEigenvalue(const SparseMT& A, const SparseMT& B, const SparseMT& C,
                                 SparseMT& /*X*/, SparseMT& /*Y*/, VT& eigvec,
                                 bool /*updateOnly*/ = false,
                                 int  /*num_constraints*/ = 0)
    {
        NT result = solveSparseQEP(A, B, C, eigvec);
        
        if (result <= NT(0) || !std::isfinite(result)) {
            return std::numeric_limits<NT>::infinity();
        }
        
        return result;
    }

    /// Find minimum positive eigenvalue for linear generalized eigenvalue problem (Sparse Version)
    /// Solves (C + t*B)v = 0 to find smallest positive t
    /// This is a special case where we have a generalized eigenvalue problem B*v = lambda*(-C)*v
    /// \param[in] C Constant coefficient matrix (sparse)
    /// \param[in] B Linear coefficient matrix (sparse)
    /// \param[out] eigvec Eigenvector corresponding to minimum positive eigenvalue
    /// \return Minimum positive t, or infinity if none exists
    NT minPosLinearEigenvalue(const SparseMT& C, const SparseMT& B, VT& eigvec)
    {
        // Solve the generalized eigenvalue problem: B*v = lambda*(-C)*v
        // We want the smallest positive eigenvalue (lambda > 0)
        std::pair<NT, NT> result = symGeneralizedProblem(C, B);
        
        NT min_pos = result.first;
        
        if (min_pos <= NT(0) || !std::isfinite(min_pos) || min_pos >= LARGE_VAL()) {
            return std::numeric_limits<NT>::infinity();
        }
        
        return min_pos;
    }

    static std::pair<NT, NT> symGeneralizedProblem(const SparseMT& A, const SparseMT& B) {
    const int n = A.rows();
    const NT tolerance = std::max(NT(1e-12),
                                   std::max(A.norm(), B.norm()) * eps() * NT(100000));
    try {
        // Convert sparse to dense for eigenvalue computation
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> dense_A(A);
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> dense_B(B);
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> neg_dense_A = -dense_A;
        
        // Solve generalized eigenvalue problem: B*v = lambda*(-A)*v
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>>
            solver(dense_B, neg_dense_A, Eigen::EigenvaluesOnly);
        
        if (solver.info() != Eigen::Success) {
            return {LARGE_VAL(), -LARGE_VAL()};
        }
        
        auto eigenvals = solver.eigenvalues();
        NT min_pos = LARGE_VAL();
        NT max_neg = -LARGE_VAL();
        
        // Binary search for the transition from negative to positive eigenvalues
        // Find the rightmost negative eigenvalue (max_neg)
        int left = 0, right = n - 1;
        int last_negative = -1;
        
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (eigenvals(mid) < -tolerance) {
                last_negative = mid;
                left = mid + 1;  // Search right for larger negative values
            } else {
                right = mid - 1;
            }
        }
        
        if (last_negative >= 0) {
            max_neg = eigenvals(last_negative);
        }
        
        // Binary search for the leftmost positive eigenvalue (min_pos)
        left = 0;
        right = n - 1;
        int first_positive = -1;
        
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (eigenvals(mid) > tolerance) {
                first_positive = mid;
                right = mid - 1;  // Search left for smaller positive values
            } else {
                left = mid + 1;
            }
        }
        
        if (first_positive >= 0) {
            min_pos = eigenvals(first_positive);
        }
        
        return {min_pos, max_neg};
    } catch (...) {
        return {LARGE_VAL(), -LARGE_VAL()};
    }
}

    // /// Check if symmetric sparse matrix M is negative definite (Sparse Version)
    // /// and return an estimate of the largest eigenvalue (most negative)
    // /// Uses sparse Cholesky factorization for efficiency
    // /// \param[in] M Symmetric sparse matrix
    // /// \return Estimated largest eigenvalue if M is negative definite, -1 or infinity otherwise
    static NT findSymEigenvalue(const SparseMT& M) {
        const int n = M.rows();
        
        // Early exit: check diagonal (sparse-friendly check)
        for (int k = 0; k < n; ++k) {
            NT diag_val = M.coeff(k, k);
            if (diag_val >= NT(0)) {
                // Return -1 to indicate not negative definite
                return NT(-1);
            }
        }
        
        // Use sparse factorization on -M
        SparseMT neg_M = -M;
        
        // Try SimplicialLLT first (faster, no pivoting)
        Eigen::SimplicialLLT<SparseMT> llt;
        llt.compute(neg_M);
        if (llt.info() == Eigen::Success) {
            // Successfully factored, -M is PD, so M is ND
            // For SimplicialLLT: -M = P' * L * L' * P
            // Get diagonal of L
            SparseMT L = llt.matrixL();
            
            // Extract diagonal values efficiently
            NT min_diag_sq = std::numeric_limits<NT>::max();
            for (int k = 0; k < L.outerSize(); ++k) {
                for (typename SparseMT::InnerIterator it(L, k); it; ++it) {
                    if (it.row() == it.col()) {
                        NT diag_val = it.value();
                        NT diag_sq = diag_val * diag_val;
                        if (diag_sq < min_diag_sq) {
                            min_diag_sq = diag_sq;
                        }
                    }
                }
            }
            return min_diag_sq;
        }
        
        // Fallback to SimplicialLDLT (handles indefinite/ill-conditioned)
        Eigen::SimplicialLDLT<SparseMT> ldlt;
        ldlt.compute(neg_M);
        if (ldlt.info() != Eigen::Success) {
            return std::numeric_limits<NT>::infinity();
        }
        
        // Extract D diagonal vector and return minimum
        return ldlt.vectorD().minCoeff();
    }
};

#endif // VOLESTI_SPARSE_EIGENVALUES_PROBLEMS_H