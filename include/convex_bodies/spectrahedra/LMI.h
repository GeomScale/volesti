// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_LMI_H
#define VOLESTI_LMI_H

#include "matrix_operations/EigenvaluesProblems.h"


template <typename MT>
struct evaluate_lmi {
    
};

template <typename NT>
struct evaluate_lmi<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> > {
public:
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    MT vectorMatrix;

    int _m, _d;

    /// Create the vectorMatrix, which has at each column the distinct elements of each A_i, i=1,...,d
    void setVectorMatrix(int const& m, int const& d, std::vector<MT> &matrices) 
    {
        _m = m;
        _d = d;
        int newM = m * (m + 1) / 2;

        // allocate memory
        vectorMatrix.setZero(newM, d);
        //a.setZero(m);

        // initialze iterator and skip A_0
        typename std::vector<MT>::iterator iter = matrices.begin();
        iter++;

        // copy elements
        int atMatrix = 0;

        for (; iter != matrices.end(); iter++, atMatrix++) {
            int i = 0;

            for (int at_row = 0; at_row < m; at_row++)
                for (int at_col = at_row; at_col < m; at_col++) {
                    vectorMatrix(i++, atMatrix) = (*iter)(at_row, at_col);
                }

        }
    }

    /// Compute  \[x_1*A_1 + ... + x_n A_n]
    /// \param[in] x Input vector
    /// \param[out] res Output matrix
    void evaluateWithoutA0(const VT &x, MT& res, bool complete_mat = false)  const {
        //#define EVALUATE_WITHOUT_A0_NAIVE
        #if defined(EVALUATE_WITHOUT_A0_NAIVE)
            res.setZero(_m, _m);
            typename std::vector<MT>::iterator it;

            int i = 0;
            it = matrices.begin();
            ++it; // skip A0
            for (; it != matrices.end(); it++, i++)
                res.noalias() += x(i) * (*it);
        #else

            VT a = vectorMatrix * x;

            double *data = res.data();
            double *v = a.data();

            int at = 0;

            // copy lower triangular
            for (int at_col = 0; at_col < _m; at_col++) {
                int col_offset = at_col * _m;
                double *target = data + col_offset + at_col;

                for (int at_row = at_col; at_row < _m; at_row++) {
                    *(target++) = *(v++);
                }
            }

            if(complete_mat) {
                v = a.data();

                // copy upper triangular
                for (int at_row = 0; at_row < _m; at_row++) {
                    double *target = data + at_row + at_row * _m;

                    for (int at_col = at_row; at_col < _m; at_col++) {
                        *target = *(v++);
                        target = target + _m;
                    }
                }
            }
        #endif
        //return true;
    }

};


/// This class handles a linear matrix inequality of the form \[A_0 +  \sum x_i A_i\]
/// A template specialization for dense Eigen matrices and vectors
/// @tparam NT Numeric Type
/// @tparam MT Matrix Type
/// @tparam VT Vector Type
template<typename NT, typename MT, typename VT>
class LMI {
    public:

    evaluate_lmi<MT> lmi_evaluator;

    /// The matrices A_0, A_i
    std::vector<MT> matrices;

    /// The dimension of the vector x
    unsigned int d;

    /// The size of the matrices A_i
    unsigned int m;

    /// At each column keep the m*(m+1)/2 distinct elements of each matrix A_i, i=1,...,d
    MT vectorMatrix;

    LMI(){}

    /// Creates A LMI object
    /// \param[in] matrices The matrices A_0, A_i
    LMI(std::vector<MT>& matrices) {
        typename std::vector<MT>::iterator it = matrices.begin();

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
    std::vector<MT> getMatrices() const {
        return matrices;
    }

    /// \returns The size of the matrices
    unsigned int sizeOfMatrices() const {
        return m;
    }

    /// Evaluate A_0 + \[A_0 + \sum x_i A_i \]
    /// \param[in] x The input vector
    /// \param[out] ret The output matrix
    void evaluate(VT const & x, MT& ret, bool complete_mat = false) const {
        lmi_evaluator.evaluateWithoutA0(x, ret, complete_mat);

        // add A0
        ret += matrices[0];
    }

    /// Compute  \[x_1*A_1 + ... + x_n A_n]
    /// \param[in] x Input vector
    /// \param[out] res Output matrix
    void evaluateWithoutA0(const VT& x, MT& res, bool complete_mat = false)  const {
        lmi_evaluator.evaluateWithoutA0(x, res, complete_mat);
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: lmi(p)*e = 0, e != 0
    /// \param[out] ret The normalized gradient of the determinant of the LMI at p
    void normalizedDeterminantGradient(VT r, VT const& e, VT &ret) const {
        NT* ret_data = ret.data();
        NT sum_sqqrt_sq = NT(0);
        for (int i = 0; i < d; i++) {
            // todo, use iterators
            *ret_data = e.dot(matrices[i+1].template selfadjointView< Eigen::Lower >() * e);
            sum_sqqrt_sq += (*ret_data) * (*ret_data);
            ret_data++;
        }

        //normalize
        ret /= std::sqrt(sum_sqqrt_sq);
    }

    /// \param i An indicator to a matrix
    /// \return Pointer to A_i
    MT* const getMatrix(const int i) {
        return &(matrices[i]);
    }

    MT get_A0() {
        return matrices[0];
    }

    void set_A0(MT const& A0) {
        matrices[0] = A0;
    }

    /// Prints the matrices A0, ..., An
    void print() const {
        int i = 0;

        for (auto iter = matrices.begin(); iter != matrices.end(); iter++, i++) {
            std::cout << "A" << i << "\n";
            std::cout << *iter << "\n\n";
        }
    }

    /// check if the matrix is negative semidefinite
    /// \param matrix a matrix
    /// \return Pointer to A_i
    bool isNegativeSemidefinite(MT const & matrix ) const {
        EigenvaluesProblems<NT, MT, VT> eigs;
        NT eival = eigs.findSymEigenvalue(matrix);
        return eival <= 0;
    }

    /// evaluate LMI(pos) and check if its negative semidefinite
    /// \param pos a vector of our current position
    /// \return true is LMI(pos) is negative semidefinite
    bool isNegativeSemidefinite(VT const & pos) const {
        MT mat;
        mat.setZero(m, m);
        
        evaluate(pos, mat, true);
        return isNegativeSemidefinite(mat);
    }

};

#endif //VOLESTI_LMI_H
