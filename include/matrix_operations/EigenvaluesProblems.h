// VolEsti (volume computation and sampling library)

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and modified by Huu Phuoc Le as part of Google Summer of Code 2022 program
// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_EIGENVALUESPROBLEMS_H
#define VOLESTI_EIGENVALUESPROBLEMS_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Spectra library for eigenvalue problems
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/DenseGenMatProd.h>

/// Solver for various eigenvalue problems arising in convex optimization
/// Provides methods for quadratic eigenvalue problems (QEP), generalized eigenvalue
/// problems, and negative definiteness checks via eigenvalue analysis.
/// \tparam NT Numeric type for scalar values
/// \tparam MT Matrix type (defaults to Eigen::Matrix with dynamic dimensions)
/// \tparam VT Vector type (defaults to Eigen::Matrix column vector with dynamic size)
template<typename NT,
         typename MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>,
         typename VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>>
class EigenvaluesProblems {
private:
    /// Machine epsilon for numeric type NT
    static constexpr NT eps() { return std::numeric_limits<NT>::epsilon(); }
    
    /// Square root of machine epsilon (used for moderate tolerance checks)
    static constexpr NT sqrt_eps() { return std::sqrt(eps()); }
    
    /// Largest representable value for numeric type NT
    static constexpr NT LARGE_VAL() { return std::numeric_limits<NT>::max(); }
    
    /// Machine epsilon alias for consistency with existing code
    static constexpr NT EPS() { return std::numeric_limits<NT>::epsilon(); }

    /// Operator for structured quadratic eigenvalue problem (QEP) linearization
    /// Implements the matrix-vector product for the linearized system C1^{-1}C0
    /// where the original QEP is: (A + λB + λ²C)v = 0
    class QEPOperator {
    private:
        const MT& A_;  ///< Coefficient matrix A in QEP
        const MT& B_;  ///< Coefficient matrix B in QEP
        const MT& C_;  ///< Coefficient matrix C in QEP (must be positive definite)
        const int n_;  ///< Dimension of original problem
        mutable Eigen::LLT<MT> lltC_;    ///< Cholesky decomposition of C (LLT)
        mutable Eigen::LDLT<MT> ldltC_;  ///< LDLT decomposition of C (fallback)
        mutable bool chol_computed_ = false;  ///< Flag indicating if decomposition is computed
        mutable bool use_llt_ = false;        ///< Flag indicating which decomposition to use
        
    public:
        using Scalar = NT;  ///< Scalar type required by Spectra library

        /// Construct QEP operator from coefficient matrices
        /// \param[in] A Constant term matrix
        /// \param[in] B Linear term matrix
        /// \param[in] C Quadratic term matrix (should be positive definite)
        QEPOperator(const MT& A, const MT& B, const MT& C)
            : A_(A), B_(B), C_(C), n_(A.rows()) {}

        /// Number of rows in linearized system (twice original dimension)
        int rows() const { return 2 * n_; }
        
        /// Number of columns in linearized system (twice original dimension)
        int cols() const { return 2 * n_; }

        /// Perform matrix-vector product y = Op * x for Spectra eigenvalue solver
        /// Implements the linearized QEP system operator using efficient factorizations
        /// \param[in] x_in Input vector of size 2*n
        /// \param[out] y_out Output vector of size 2*n
        void perform_op(const NT* x_in, NT* y_out) const {
            Eigen::Map<const VT> x(x_in, 2*n_);
            Eigen::Map<VT> y(y_out, 2*n_);
            
            // Compute factorization once on first call (lazy initialization)
            // Try LLT first (faster for positive definite matrices)
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
            rhs.noalias() = B_ * x.head(n_);
            rhs += x.tail(n_);
            
            // Solve C * y.head(n) = -rhs using precomputed factorization
            if (use_llt_) {
                y.head(n_).noalias() = lltC_.solve(-rhs);
            } else if (ldltC_.info() == Eigen::Success) {
                y.head(n_).noalias() = ldltC_.solve(-rhs);
            } else {
                // Fallback to QR decomposition (extremely rare, for robustness)
                y.head(n_).noalias() = C_.colPivHouseholderQr().solve(-rhs);
            }
            
            // Set bottom half of output: y.tail(n) = x.head(n)
            y.tail(n_) = x.head(n_);
        }
    };

    /// Solve quadratic eigenvalue problem to find smallest positive parameter
    /// Solves (B0 + t*B1 + t²*B2)v = 0 for the smallest positive t
    /// Uses Spectra library with linearization approach via QEPOperator
    /// \param[in] B0 Constant coefficient matrix
    /// \param[in] B1 Linear coefficient matrix
    /// \param[in] B2 Quadratic coefficient matrix
    /// \param[out] eigvec Eigenvector corresponding to smallest positive t
    /// \return Smallest positive value of t, or infinity if none exists
    static NT solveQEP(const MT& B0, const MT& B1, const MT& B2, VT& eigvec) {
    const int m = B0.rows();
    if (m == 0) return std::numeric_limits<NT>::infinity();
    
    try {
        QEPOperator op(B0, B1, B2);
        
        // Empirical choice: request 2 eigenvalues to improve chances of finding valid positive t
        // for faster convergence use nev = 1, but more unstable, may miss valid eigenvalue
        const int nev = 1; 
        const int ncv = std::min(10, 2*m);
        
        Spectra::GenEigsSolver<QEPOperator> solver(op, nev, ncv);
        solver.init();
        
        int nconv = solver.compute(Spectra::SortRule::LargestMagn, 1500, NT(1e-10));
        
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
    /// Find minimum positive eigenvalue for quadratic eigenvalue problem
    /// Solves (A + t*B + t²*C)v = 0 to find smallest positive t
    /// Used in convex optimization for barrier function step size computation
    /// \param[in] A Constant coefficient matrix
    /// \param[in] B Linear coefficient matrix
    /// \param[in] C Quadratic coefficient matrix
    /// \param[in,out] X Auxiliary matrix (unused, kept for API compatibility)
    /// \param[in,out] Y Auxiliary matrix (unused, kept for API compatibility)
    /// \param[out] eigvec Eigenvector corresponding to minimum positive eigenvalue
    /// \param[in] updateOnly Unused flag (kept for API compatibility)
    /// \param[in] num_constraints Unused parameter (kept for API compatibility)
    /// \return Minimum positive t satisfying the QEP, or infinity if none exists
    NT minPosQuadraticEigenvalue(const MT& A, const MT& B, const MT& C,
                                 MT& /*X*/, MT& /*Y*/, VT& eigvec,
                                 bool /*updateOnly*/ = false,
                                 int  /*num_constraints*/ = 0)
    {
        
        NT result = solveQEP(A, B, C, eigvec);
        
        if (result <= NT(0) || !std::isfinite(result)) {
            return std::numeric_limits<NT>::infinity();
        }
        
        return result;
    }

    /// Solve generalized symmetric eigenvalue problem B*v = λ*(-A)*v
    /// Finds both minimum positive and maximum negative eigenvalues simultaneously
    /// \param[in] A Symmetric matrix (multiplied by -1 in problem formulation)
    /// \param[in] B Symmetric matrix
    /// \return Pair (min_positive_eigenvalue, max_negative_eigenvalue)
    static std::pair<NT, NT> symGeneralizedProblem(const MT& A, const MT& B) {
        if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
            return {LARGE_VAL(), -LARGE_VAL()};
        }

        // retrieve corresponding eigenvector
        eigenvector = ges.eigenvectors().col(index);
#elif defined(SPECTRA_EIGENVALUES_SOLVER)
        // Transform the problem to a standard eigenvalue problem and use the general eigenvalue solver of Spectra

        // This makes the transformation to standard eigenvalue problem. See class for more info.
        // We have the generalized problem  A + lB, or Av = -lBv
        // This class computes the matrix product vector Mv, where M = -B * A^[-1]

        DenseProductMatrix<NT> M(&B, &A,true);

        // This parameter is for Spectra. It must be larger than #(requested eigenvalues) + 2
        // and smaller than the size of matrix;
        int ncv = 3;

        // Prepare to solve Mx = (1/l)x
        // we want the smallest positive eigenvalue in the original problem,
        // so in this the largest positive eigenvalue;
        Spectra::GenEigsSolver<NT, Spectra::LARGEST_REAL, DenseProductMatrix<NT> > eigs(&M, 1, ncv);

        // compute
        eigs.init();
        eigs.compute();

        //retrieve result and invert to get required eigenvalue of the original problem
        if (eigs.info() != Spectra::SUCCESSFUL) {
            eigenvector.setZero(A.rows());
            return NT(0);
        }

        lambdaMinPositive = 1/((eigs.eigenvalues())(0).real());

        // retrieve corresponding eigenvector
        int matrixDim = A.rows();
        eigenvector.resize(matrixDim);
        for (int i = 0; i < matrixDim; i++)
            eigenvector(i) =  (eigs.eigenvectors()).col(0)(i);

#elif defined(ARPACK_EIGENVALUES_SOLVER)
        // Transform the problem to a standard eigenvalue problem and use the general eigenvalue solver of ARPACK++

        // This makes the transformation to standard eigenvalue problem. See class for more info.
        // We have the generalized problem  A + lB, or Av = -lBv
        // This class computes the matrix product vector Mv, where M = -B * A^[-1]

        DenseProductMatrix<NT> M(&B, &A,true);

        // Creating an eigenvalue problem and defining what we need:
        // the  eigenvector of A with largest real.
        ARNonSymStdEig<NT, DenseProductMatrix<NT> >

        dprob(A.cols(), 1, &M, &DenseProductMatrix<NT>::MultMv, std::string ("LR"), 8<A.rows() ? 8 : A.rows(), 0.000);//, 100*3);

        // compute
        if (dprob.FindEigenvectors() == 0) {
            std::cout << "Failed\n";
            // if failed with default (and fast) parameters, try with stable (and slow)
            dprob.ChangeNcv(A.cols()/10);
            if (dprob.FindEigenvectors() == 0) {
                std::cout << "\tFailed Again\n";
                return NT(0);
        
        const int n = A.rows();
        const NT tolerance = std::max(NT(1e-12), 
                                    std::max(A.norm(), B.norm()) * eps() * NT(100000));
        
        try {
            MT neg_A = -A;
            Eigen::GeneralizedSelfAdjointEigenSolver<MT> solver(B, neg_A, Eigen::EigenvaluesOnly);
            
            if (solver.info() != Eigen::Success) {
                return {LARGE_VAL(), -LARGE_VAL()};
            }
            
            auto eigenvals = solver.eigenvalues();
            
            NT min_pos = LARGE_VAL();
            NT max_neg = -LARGE_VAL();
            
            // Single forward pass (eigenvalues sorted increasing)
            // Keep updating max_neg, then grab first min_pos and exit
            for (int i = 0; i < n; ++i) {
                NT lambda = eigenvals(i);
                
                if (std::abs(lambda) < tolerance) continue;
                
                if (lambda < -tolerance) {
                    max_neg = lambda;  // Keep updating (moving toward zero)
                } else if (lambda > tolerance) {
                    min_pos = lambda;  // First positive = minimum positive
                    break;  // Done! Have both values
                }
            }
            
            return {min_pos, max_neg};
        } catch (...) {
            return {LARGE_VAL(), -LARGE_VAL()};
        }
    }

    /// Check if symmetric matrix M is negative definite
    /// and return an estimate of the largest eigenvalue (most negative)
    /// \param[in] M Symmetric matrix
    /// \return Estimated largest eigenvalue if M is negative definite, infinity otherwise
    static NT findSymEigenvalue(const MT& M) {
        const int n = M.rows();
        
        // Early exit: check diagonal
        for (int i = 0; i < n; ++i) {
            if (M.coeff(i, i) >= NT(0)) {
                // Return -1 to indicate not negative definite
                return -1;
            }
        }
        
        // Try Cholesky first (faster, no pivoting)
        Eigen::LLT<MT> llt;
        llt.compute(-M);
        if (llt.info() == Eigen::Success) {
            // Successfully factored, -M is PD, so M is ND
            // Get diagonal elements from L: L is lower triangular
            // For LLT: -M = L*L^T, diagonal of D would be L.diagonal()^2
            Eigen::VectorXd diag = llt.matrixLLT().diagonal();
            // Return min of squared diagonal (approximates eigenvalue bound)
            return diag.array().square().minCoeff();
        }
        
        // Fallback to LDLT (handles indefinite/ill-conditioned)
        Eigen::LDLT<MT> ldlt;
        ldlt.compute(-M);
        if (ldlt.info() != Eigen::Success) {
            return std::numeric_limits<NT>::infinity();
        }
        return ldlt.vectorD().minCoeff();
    }
};

#endif // VOLESTI_EIGENVALUESPROBLEMS_H