//
// Created by panagiotis on 2/22/20.
//

#ifndef VOLESTI_EIGENDENSEMATRIX_H
#define VOLESTI_EIGENDENSEMATRIX_H

/// A wrap class to use Eigen dense matrices when solving Eigenvalue problems with ARPACK++
/// \tparam NT Numeric Type
template<class NT>
class EigenDenseMatrix {
public:

    /// Eigen matrix type
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    /// Eigen vector type
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    /// The matrix
    MT const * M;

    /// number of columns
    int n;
    /// number of rows
    int m;

    /// \return Number of rows
    int nrows() { return m;}

    /// \return Number of columns
    int ncols() { return n;}

    /// \return Number of rows
    int rows() { return m;}

    /// \return Number of columns
    int cols() { return n;}

    /// Required by ARPACK++ : Multiplies the matrix with vector v
    /// \param[in] v The input vector, for example double*
    /// \param[out] w The result of M*v
    void MultMv(NT* v, NT* w) {
        // Declaring the vectors like this, we don't copy the values of v and after to w
        Eigen::Map<VT> _v(v, m);
        Eigen::Map<VT> _w(w, m);

        _w = *M * _v;
    }

    /// Required by ARPACK++ : Multiplies the matrix with vector v
    /// \param[in] v The input vector, for example double*
    /// \param[out] w The result of M*v
    void perform_op(NT* v, NT* w) {
        // Declaring the vectors like this, we don't copy the values of v and after to w
        Eigen::Map<VT> _v(v, m);
        Eigen::Map<VT> _w(w, m);

        _w = *M * _v;
    }


    /// Constructs an object
    /// \param[in] M An Eigen Matrix
    EigenDenseMatrix(MT const * M) {
        this->M = M;
        n = M->cols();
        m = M->rows();
    }

};
#endif //VOLESTI_EIGENDENSEMATRIX_H
