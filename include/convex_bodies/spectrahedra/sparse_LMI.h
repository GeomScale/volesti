// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SPARSE_LMI_H
#define VOLESTI_SPARSE_LMI_H

#include "matrix_operations/SparseEigenvaluesProblems.h"
#include <Eigen/Sparse>


template <typename MT>
struct evaluate_sparse_lmi {
    
};

template <typename NT>
struct evaluate_sparse_lmi<Eigen::SparseMatrix<NT> > {
public:
    typedef Eigen::SparseMatrix<NT> SparseMT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DenseMT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    DenseMT vectorMatrix;

    int _m, _d;

    /// Create the vectorMatrix, which has at each column the distinct elements of each A_i, i=1,...,d
    /// For sparse matrices, we iterate over stored values directly for efficiency
    void setVectorMatrix(int const& m, int const& d, std::vector<SparseMT> &matrices) 
    {
        _m = m;
        _d = d;
        int newM = m * (m + 1) / 2;

        // allocate memory for dense vectorMatrix
        vectorMatrix.setZero(newM, d);

        // initialize iterator and skip A_0
        typename std::vector<SparseMT>::iterator iter = matrices.begin();
        iter++;

        // copy elements from sparse matrices
        int atMatrix = 0;

        for (; iter != matrices.end(); iter++, atMatrix++) {
            // Iterate over stored values efficiently (O(nnz) instead of O(m²))
            for (int k = 0; k < iter->outerSize(); ++k) {
                for (typename SparseMT::InnerIterator it(*iter, k); it; ++it) {
                    int row = it.row();
                    int col = it.col();
                    NT val = it.value();
                    
                    // Map (row, col) to vectorMatrix index
                    // We store upper triangle: positions where col >= row
                    int vrow, vcol;
                    if (col >= row) {
                        vrow = row;
                        vcol = col;
                    } else {
                        // Stored in lower triangle, swap for upper triangle
                        vrow = col;
                        vcol = row;
                    }
                    
                    // Calculate linear index in upper triangle storage
                    // Position = row * m - row*(row-1)/2 + (col - row)
                    int idx = vrow * m - vrow * (vrow - 1) / 2 + (vcol - vrow);
                    
                    vectorMatrix(idx, atMatrix) = val;
                }
            }
        }
    }

    /// Compute  \[x_1*A_1 + ... + x_n A_n] for sparse matrices
    /// \param[in] x Input vector
    /// \param[out] res Output sparse matrix
    void evaluateWithoutA0(const VT &x, SparseMT& res, bool complete_mat = false)  const {
        // Use dense vectorMatrix multiplication for efficiency
        VT a = vectorMatrix * x;

        // Prepare triplets for sparse matrix construction
        std::vector<Eigen::Triplet<NT>> triplets;
        triplets.reserve(_m * (_m + 1));

        const NT* v = a.data();

        // Fill in column-major order matching dense version
        int at = 0;
        for (int at_col = 0; at_col < _m; at_col++) {
            for (int at_row = at_col; at_row < _m; at_row++) {
                NT val = v[at++];
                triplets.push_back(Eigen::Triplet<NT>(at_row, at_col, val));
                
                // If symmetric completion requested and not on diagonal
                if (complete_mat && at_row != at_col) {
                    triplets.push_back(Eigen::Triplet<NT>(at_col, at_row, val));
                }
            }
        }

        // Build sparse matrix from triplets
        res.resize(_m, _m);
        res.setFromTriplets(triplets.begin(), triplets.end());
        res.makeCompressed();
    }
};


/// This class handles a linear matrix inequality of the form \[A_0 +  \sum x_i A_i\]
/// A template specialization for sparse Eigen matrices and vectors
/// @tparam NT Numeric Type
/// @tparam SparseMT Sparse Matrix Type (Eigen::SparseMatrix)
/// @tparam VT Vector Type
template<typename NT, typename SparseMT, typename VT>
class SparseLMI {
public:

    evaluate_sparse_lmi<SparseMT> lmi_evaluator;

    /// The matrices A_0, A_i
    std::vector<SparseMT> matrices;

    /// The dimension of the vector x
    unsigned int d;

    /// The size of the matrices A_i
    unsigned int m;

    /// At each column keep the m*(m+1)/2 distinct elements of each matrix A_i, i=1,...,d
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> vectorMatrix;

    SparseLMI(){}

    /// Creates A LMI object with sparse matrices
    /// \param[in] matrices The matrices A_0, A_i (sparse)
    SparseLMI(std::vector<SparseMT>& matrices) {
        typename std::vector<SparseMT>::iterator it = matrices.begin();

        while (it!=matrices.end()) {
            this->matrices.push_back(*it);
            it++;
        }

        d = matrices.size() - 1;
        m = matrices[0].rows();

        lmi_evaluator.setVectorMatrix(m, d, matrices);
    }

    /// \returns The dimension of vector x
    unsigned int dimension() const {
        return d;
    }

    /// \return The matrices A0, A1, ..., Ad
    std::vector<SparseMT> getMatrices() const {
        return matrices;
    }

    /// \returns The size of the matrices
    unsigned int sizeOfMatrices() const {
        return m;
    }

    /// Evaluate A_0 + \[A_0 + \sum x_i A_i \]
    /// \param[in] x The input vector
    /// \param[out] ret The output sparse matrix
    void evaluate(VT const & x, SparseMT& ret, bool complete_mat = false) const {
        lmi_evaluator.evaluateWithoutA0(x, ret, complete_mat);

        // add A0 using sparse addition
        ret = ret + matrices[0];
    }

    /// Compute  \[x_1*A_1 + ... + x_n A_n]
    /// \param[in] x Input vector
    /// \param[out] res Output sparse matrix
    void evaluateWithoutA0(const VT& x, SparseMT& res, bool complete_mat = false)  const {
        lmi_evaluator.evaluateWithoutA0(x, res, complete_mat);
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] r Input parameter (kept for API compatibility)
    /// \param[in] e Input vector: lmi(p)*e = 0, e != 0
    /// \param[out] ret The normalized gradient of the determinant of the LMI at p
    void normalizedDeterminantGradient(VT r, VT const& e, VT &ret) const {
        NT* ret_data = ret.data();
        NT sum_sqrt_sq = NT(0);
        
        for (int i = 0; i < d; i++) {
            // For sparse matrices, compute e^T * A_i * e
            // Convert to dense for reliable symmetric view (matching dense version's selfadjointView)
            Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> dense_mat = matrices[i+1];
            *ret_data = e.dot(dense_mat.template selfadjointView<Eigen::Lower>() * e);
            
            sum_sqrt_sq += (*ret_data) * (*ret_data);
            ret_data++;
        }

        //normalize
        ret /= std::sqrt(sum_sqrt_sq);
    }

    /// \param i An indicator to a matrix
    /// \return Pointer to A_i
    SparseMT* const getMatrix(const int i) {
        return &(matrices[i]);
    }

    SparseMT get_A0() {
        return matrices[0];
    }

    void set_A0(SparseMT const& A0) {
        matrices[0] = A0;
    }

    /// Prints the matrices A0, ..., An
    void print() const {
        int i = 0;

        for (auto iter = matrices.begin(); iter != matrices.end(); iter++, i++) {
            std::cout << "A" << i << "\n";
            std::cout << Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>(*iter) << "\n\n";
        }
    }

    /// check if the matrix is negative definite
    /// \param matrix a sparse matrix
    /// \return true if matrix is negative definite
    bool isNegativeDefinite(SparseMT const & matrix) const {
        const NT tol = NT(1e-10) * matrix.norm();

        SparseEigenvaluesProblems<NT, SparseMT, VT> eigs;
        NT eival = eigs.findSymEigenvalue(matrix);  
        return eival >= tol;
    }

    /// evaluate LMI(pos) and check if its negative definite
    /// \param pos a vector of our current position
    /// \return true if LMI(pos) is negative definite
    bool isNegativeDefinite(VT const & pos) const {
        SparseMT mat;
        mat.resize(m, m);
        mat.setZero();
        
        evaluate(pos, mat, true);
        return isNegativeDefinite(mat);
    }

};

#endif //VOLESTI_SPARSE_LMI_H