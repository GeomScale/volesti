#ifndef VOLESTI_EIGENVALUES_CORRELATION_H
#define VOLESTI_EIGENVALUES_CORRELATION_H

//#define EIGEN_EIGENVALUES_SOLVER // Eigen
#define SPECTRA_EIGENVALUES_SOLVER // Spectra 
//#define ARPACK_EIGENVALUES_SOLVER // ARPACK++ 

#include <../../external/Spectra/include/Spectra/SymEigsSolver.h>
#include "DenseProductMatrix.h"
#include "EigenDenseMatrix.h"

#include "../../external/Spectra/include/Spectra/SymGEigsSolver.h"
#include "../../external/Spectra/include/Spectra/GenEigsSolver.h"

/// \tparam NT, MT, VT
template<typename NT, typename MT, typename VT>
class EigenvaluesCorrelation {
public:

    typedef std::pair<NT, NT> NTpair;

    #if defined(EIGEN_EIGENVALUES_SOLVER) || defined (SPECTRA_EIGENVALUES_SOLVER)
    typedef typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType CVT;
    #elif defined(ARPACK_EIGENVALUES_SOLVER)
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> CVT;
    #endif

    // Using LDLT decomposition: more numerically stable for singular matrices
    bool isPositiveSemidefinite(MT const &A) {
        Eigen::LDLT<MT> A_ldlt(A);
        if (A_ldlt.info() != Eigen::NumericalIssue && A_ldlt.isPositive())
            return true;
        return false;
    }


    /// Find the smallest eigenvalue of mat
    /// \param mat a symmetric matrix
    /// \return the smallest eigenvalue of mat
    NT smallestEigenvalue(MT const & mat) const {
#if defined(SPECTRA_EIGENVALUES_SOLVER)
        Spectra::DenseSymMatProd<NT> M(-mat);
        int ncv = M.rows()/10 + 5;
        if (ncv > M.rows()) ncv = M.rows();
        Spectra::SymEigsSolver<NT, Spectra::SmallestAlge, Spectra::DenseSymMatProd<NT>> eigs(&M, 1, ncv);
        // compute
        eigs.init();
        eigs.compute(50000);
        if(eigs.info() == Spectra::SUCCESSFUL)
            return -eigs.eigenvalues()(0);
#else
            EigenDenseMatrix<NT> M(&mat);
            Eigen::SelfAdjointEigenSolver<MT> solver;
            solver.compute(M, Eigen::EigenvaluesOnly);
            return solver.eigenvalues().minCoeff();
#endif
    }

    /// Minimum positive and maximum negative eigenvalues of the generalized eigenvalue problem A - lB
    /// \param[in] A: symmetric positive definite matrix
    /// \param[in] B: symmetric matrix
    /// \return The pair (minimum positive, maximum negative) of eigenvalues
    NTpair symGeneralizedProblem(MT const & A, MT const & B) {
        int matrixDim = A.rows();
        int ncv = 5 < matrixDim ? 5 : matrixDim; // convergence speed

        Spectra::DenseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(A);

        // Request the smallest and largest positive eigenvalues
        Spectra::SymGEigsSolver<NT, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, ncv);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        if (geigs.info() != Spectra::SUCCESSFUL)
            return {NT(0), NT(0)};

        double lambdaMinPositive, lambdaMaxNegative;

        Eigen::VectorXd evalues = geigs.eigenvalues();

        // get the eigenvalues of the original problem
        lambdaMinPositive = 1 / evalues(0);
        lambdaMaxNegative = 1 / evalues(1);

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /// Minimum positive eigenvalue of the generalized eigenvalue problem A - lB
    /// \param[in] A: symmetric positive definite matrix
    /// \param[in] B: symmetric matrix
    /// \return The minimum positive eigenvalue and the corresponding eigenvector
    NT minPosLinearEigenvalue(MT const & A, MT const & B, VT &eigvec) {
        int matrixDim = A.rows();
        int ncv = 15 < matrixDim ? 15 : matrixDim; // convergence speed

        Spectra::DenseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(A);

        // Requesting the largest generalized eigenvalues
        Spectra::SymGEigsSolver<NT, Spectra::LARGEST_ALGE,  Spectra::DenseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY> 
            geigs(&op, &Bop, 1, ncv);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        VT evalues;
        if (geigs.info() == Spectra::SUCCESSFUL) {
            evalues = geigs.eigenvalues();
            eigvec = geigs.eigenvectors().col(0);
        }

        double lambdaMinPositive = 1 / evalues(0);

        return lambdaMinPositive;
    }

    #if defined(EIGEN_EIGENVALUES_SOLVER)
        // use the Generalized eigenvalue solver of Eigen

        // compute generalized eigenvalues with Eigen solver
        Eigen::GeneralizedEigenSolver<MT> ges(A, -B);

        // retrieve minimum positive eigenvalue
        typename Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();
        int index = 0;

        for (int i = 0; i < alphas.rows(); i++) {

            if (betas(i) == 0 || alphas(i).imag() != 0)
                continue;

            double lambda = alphas(i).real() / betas(i);
            if (lambda > 0 && lambda < lambdaMinPositive) {
                lambdaMinPositive = lambda;
                index = i;
            }
        }

        // retrieve corresponding eigenvector
        eigenvector = ges.eigenvectors().col(index);
};

#endif //VOLESTI_EIGENVALUES_CORRELATION_H
