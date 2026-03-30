#ifndef FULL_DIMENSIONAL_POLYTOPE_HPP
#define FULL_DIMENSIONAL_POLYTOPE_HPP

#include <Eigen/Eigen>
#include <tuple>
#include <stdexcept>
#include <cstring>

#include "SuiteSparseQR_C.h"
#include "SuiteSparseQR.hpp"


/**
 * Compute the full dimensional polytope P' = {x | A'x <= b'} from the polytope P = {x | Ax <= b & Aeq*x = beq}
 * 
 * @tparam NT Number type (default: double)
 * @tparam SpMT Sparse matrix type (default: Eigen::SparseMatrix<NT, Eigen::ColMajor>)
 * @tparam MT Dense matrix type (default: Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>)
 * @tparam VT Vector type (default: Eigen::Matrix<NT, Eigen::Dynamic, 1>)
 * @param Aeq Sparse matrix of equality constraints (m x n)
 * @param beq Right-hand side vector for equality constraints (m x 1)  
 * @param A Dense matrix of inequality constraints (p x n)
 * @param b Right-hand side vector for inequality constraints (p x 1)
 * 
 * @return std::tuple<MT, VT, VT, MT> containing:
 *   - A_full: Transformed inequality matrix (p x (n-rank))
 *   - b_full: Transformed inequality vector (p x 1)  
 *   - shift: Translation vector (n x 1)
 *   - N: Nullspace transformation matrix (n x (n-rank))
 */
template<typename NT = double,
         typename SpMT = Eigen::SparseMatrix<NT, Eigen::ColMajor>,
         typename MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>, 
         typename VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>>
std::tuple<MT, VT, VT, MT> compute_full_dimensional_polytope(
    const SpMT& Aeq,
    const VT& beq,
    const MT& A,
    const VT& b)
{
    cholmod_common c;
    memset(&c, 0, sizeof(cholmod_common));  // Zero-initialize the structure

    // First try solving Aeq*x = beq using Eigen's sparse QR
    Eigen::SparseQR<SpMT, Eigen::COLAMDOrdering<int>> Aqr;
    Aqr.compute(Aeq);
    VT shift = Aqr.solve(beq);

    // If eigen fails try with SuiteSparse
    if ((Aeq * shift - beq).cwiseAbs().maxCoeff() > 1e-09)
    {
        // Eigen failed to solve Aeq*x = beq
        // Allocate the Aeq
        cholmod_triplet *T = cholmod_allocate_triplet(Aeq.rows(), Aeq.cols(), Aeq.nonZeros(), 0, CHOLMOD_REAL, &c);
        int i = 0;
        for (int k=0; k<Aeq.outerSize(); ++k)
        {
            for (typename SpMT::InnerIterator it(Aeq,k); it; ++it)
            {
                ((int*)T->i)[i] = it.row();
                ((int*)T->j)[i] = it.col();
                ((double*)T->x)[i] = it.value();
                i++;
            }
        }
        T->nnz = Aeq.nonZeros();

        // Convert triplet form to sparse matrix
        cholmod_sparse *Aeq_chol = cholmod_triplet_to_sparse(T, Aeq.nonZeros(), &c);  
        // Create right-hand side vector beq for the CHOLMOD
        cholmod_dense *b_chol = cholmod_allocate_dense(Aeq.rows(), 1, Aeq.rows(), CHOLMOD_REAL, &c);
        double *b_values = (double *)b_chol->x;
        for (int j = 0; j < beq.size(); j++)
        {
            b_values[j] = beq.coeff(j);
        }
        // Solve the linear system Aeq * x = beq using SuiteSparseQR_C_backslash
        cholmod_dense *x = SuiteSparseQR_C_backslash(0, -4.0, Aeq_chol, b_chol, &c);

        if (!x) {
            cholmod_free_dense(&x, &c);
            cholmod_free_dense(&b_chol, &c);
            cholmod_free_triplet(&T, &c);
            cholmod_free_sparse(&Aeq_chol, &c);
            cholmod_finish(&c);
            throw std::runtime_error("Failed to solve the linear system.");
        }

        // copy the solution x
        shift.resize(Aeq.cols());
        NT *shift_data = shift.data();
        for (int j = 0; j < x->nrow; ++j) {
            *shift_data = static_cast<NT>(((double*)x->x)[j]);
            shift_data++;
        }

        cholmod_free_dense(&x, &c);
        cholmod_free_dense(&b_chol, &c);
        cholmod_free_triplet(&T, &c);
        cholmod_free_sparse(&Aeq_chol, &c);
    }

    // Shift the polytope so that Aeq * x = 0 to hold for the feasible points
    VT b_full = b - A * shift;

    // Allocate the Aeq.transpose()
    cholmod_triplet *T = cholmod_allocate_triplet(Aeq.cols(), Aeq.rows(), Aeq.nonZeros(), 0, CHOLMOD_REAL, &c);
    int i = 0;
    for (int k=0; k<Aeq.outerSize(); ++k)
    {
        for (typename SpMT::InnerIterator it(Aeq,k); it; ++it)
        {
            ((int*)T->i)[i] = it.col();
            ((int*)T->j)[i] = it.row();
            ((double*)T->x)[i] = it.value();
            i++;
        }
    }
    T->nnz = Aeq.nonZeros();  // Set the number of non-zero elements

    // Convert triplet form to sparse matrix
    cholmod_sparse *Aeq_tr = cholmod_triplet_to_sparse(T, Aeq.nonZeros(), &c);
    int ordering = 0;  // No permutation, natural order
    int64_t *E = NULL;
    int64_t econ = 0;         // Use the rank of Aeq_tr
    double tol = -4.0;        // default tolerance

    // Call SuiteSparseQR_C_factorize
    SuiteSparseQR_C_factorization *QR = SuiteSparseQR_C_factorize(ordering, tol, Aeq_tr, &c);

    if (!QR) {
        cholmod_free_sparse(&Aeq_tr, &c);
        cholmod_free_triplet(&T, &c);
        SuiteSparseQR_C_free(&QR, &c);
        cholmod_finish(&c);
        throw std::runtime_error("QR factorization failed.");
    }

    int rank = c.SPQR_istat[4]; // rank of Aeq.transpose()

    // Form the identity matrix
    cholmod_dense *I = cholmod_zeros(Aeq_tr->nrow, Aeq_tr->nrow, CHOLMOD_REAL, &c);
    for (int j = 0; j < Aeq_tr->nrow; j++) {
        ((double*)I->x)[j * Aeq_tr->nrow + j] = 1.0;  // Set diagonal elements to 1
    }

    // Output matrix to hold Q
    cholmod_dense *Q = NULL;

    // Retrieve Q by multiplying Q with the identity matrix I
    Q = SuiteSparseQR_C_qmult(1, // 1 for Q * I (to get Q)
                                QR, I, &c);
    if (!Q) {
        cholmod_free_sparse(&Aeq_tr, &c);
        cholmod_free_dense(&Q, &c);
        cholmod_free_dense(&I, &c);
        cholmod_free_triplet(&T, &c);
        SuiteSparseQR_C_free(&QR, &c);
        cholmod_finish(&c);
        throw std::runtime_error("Failed to retrieve Q matrix.");
    }
    int Q_rows = Q->nrow;
    int Q_cols = Q->ncol;
    double* Q_data = (double*)Q->x;

    // Map the data from CHOLMOD's Q to Eigen Matrix
    Eigen::Map<MT> QQ(Q_data, Q_rows, Q_cols);
    // Take the last n-r columns of Q to derive the right nullspace of Aeq.transpose()
    MT N = QQ.block(0, rank, QQ.rows(), QQ.cols() - rank);
    // Project the polytope to the null space
    MT A_full = A * N;

    cholmod_free_sparse(&Aeq_tr, &c);
    cholmod_free_dense(&Q, &c);
    cholmod_free_dense(&I, &c);
    cholmod_free_triplet(&T, &c);
    SuiteSparseQR_C_free(&QR, &c);
    cholmod_finish(&c);

    return std::make_tuple(A_full, b_full, shift, N);
}

#endif // FULL_DIMENSIONAL_POLYTOPE_HPP
