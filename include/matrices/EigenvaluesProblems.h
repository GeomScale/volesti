//
// Created by panagiotis on 2/22/20.
//

#ifndef VOLESTI_EIGENVALUESPROBLEMS_H
#define VOLESTI_EIGENVALUESPROBLEMS_H

#include "Spectra/MatOp/DenseCholesky.h"
#include "Spectra/MatOp/DenseSymMatProd.h"
#include "Spectra/SymGEigsSolver.h"

/// Solve eigenvalues problems
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class EigenvaluesProblems {

};


/// A specialization of the template class EigenvaluesProblems for dense Eigen matrices and vectors.
/// \tparam NT
template<typename NT>
class EigenvaluesProblems<NT, Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<NT,Eigen::Dynamic,1> > {
public:
    /// The type for Eigen Matrix
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// The type of a pair of NT
    typedef std::pair<NT, NT> NTpair;


    /// Find the minimum positive and maximum negative eigenvalues of the generalized eigenvalue
    /// problem A + lB, where A, B symmetric and A negative definite.
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \return The pair (minimum positive, maximum negative) of eigenvalues
    NTpair symGeneralizedProblem(MT const & A, MT const & B) {

        int matrixDim = A.rows();

        // Spectra solves Xv=lYv, where Y positive definite
        // Set X = B, Y=-A. Then, the eigenvalues we want are the minimum negative
        // and maximum positive eigenvalues of Xv=lYv.

        // Construct matrix operation object using the wrapper classes provided by Spectra
        Spectra::DenseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(-A);

        // Construct generalized eigen solver object
        // requesting the minmum negative and largest positive eigenvalues
        Spectra::SymGEigsSolver<NT, Spectra::BOTH_ENDS, Spectra::DenseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 2, 5 < matrixDim ? 5 : matrixDim);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        if (geigs.info() != Spectra::SUCCESSFUL)
            return {NT(0), NT(0)};

        Eigen::VectorXd evalues;
        double lambdaMinPositive, lambdaMaxNegative;

        evalues = geigs.eigenvalues();

        // get the eigenvalues of the original problem
        lambdaMinPositive = 1 / evalues(0);
        lambdaMaxNegative = 1 / evalues(1);

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /// Find the minimum positive eigenvalue of the generalized eigenvalue
    /// problem A + lB, where A, B symmetric and A negative definite.
    /// \param[in] A Input matrix
    /// \param[in] B Input matrix
    /// \param[out] eigenvector The eigenvector corresponding to the minimum positive eigenvalue
    /// \return The pair (minimum positive, maximum negative) of eigenvalues
    NT minPosSymGeneralizedProblem(MT const & A, MT const & B, VT & eigenvector) {

        int matrixDim = A.rows();

        // Spectra solves Xv=lYv, where Y positive definite
        // Set X = B, Y=-A. Then, the eigenvalue we want
        // is the maximum positive eigenvalue of Xv=lYv.

        // Construct matrix operation object using the wrapper classes provided by Spectra
        Spectra::DenseSymMatProd<NT> op(B);
        Spectra::DenseCholesky<NT> Bop(-A);

        // Construct generalized eigen solver object
        // requesting the largest positive eigenvalue
        Spectra::SymGEigsSolver<NT, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<NT>, Spectra::DenseCholesky<NT>, Spectra::GEIGS_CHOLESKY>
                geigs(&op, &Bop, 1, 5 < matrixDim ? 5 : matrixDim);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        if (geigs.info() != Spectra::SUCCESSFUL)
            return NT(0);

        Eigen::VectorXd evalues;
        double lambdaMinPositive;

        evalues = geigs.eigenvalues();

        // get the eigenvalues of the original problem
        lambdaMinPositive = 1 / evalues(0);

        eigenvector = geigs.eigenvectors().col(0);

        return lambdaMinPositive;
    }
};

#endif //VOLESTI_EIGENVALUESPROBLEMS_H
