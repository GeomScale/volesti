// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_DENSEPRODUCTMATRIX_H
#define VOLESTI_DENSEPRODUCTMATRIX_H

#define PARTIAL_LU_DECOMPOSITION

/// A wrapper class for dense Eigen matrices in Spectra and ARPACK++
/// This class will be the wrapper to use the Spectra nonsymemmetric standard eigenvalue Cx = lx solver to
/// solve a generalized eigenvalue Ax = lBx.
/// In particular, this class represents the product @f[ C = B^-1 A @f]
///
/// \tparam NT Numeric Type
template<typename NT>
class DenseProductMatrix {
public:
    /// Eigen matrix type
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    /// Eigen vector type
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// The number of rows
    int _rows;
    /// The number of cols
    int _cols;

    /// Pointer to matrix A
    MT const *A;
    /// Pointer to matrix B
    MT const *B;

    /// The decomposition we will use
    /// If PARTIAL_LU_DECOMPOSITION is defined, use the Eigen partial LU decomposition,
    /// otherwise full. The partial is faster but assumes that the matrix has full rank.
#if defined(PARTIAL_LU_DECOMPOSITION)
    typedef Eigen::PartialPivLU<MT> Decomposition;
#else
    typedef Eigen::FullPivLU<MT> Decomposition;
#endif

    /// The LU decomposition of B
    Decomposition Blu;

    /// Constructs an object of this class and computes the LU decomposition of B.
    ///
    /// \param[in] A The matrix A
    /// \param[in] B The matrix B
    DenseProductMatrix(MT const *A, MT const *B) : A(A), B(B) {
        Blu = Decomposition(*B);
        _rows = A->rows();
        _cols = B->cols();
    }

    ///Required by Spectra
    /// \return The number of rows
    int rows() {
        return _rows;
    }

    ///Required by Spectra
    /// \return The number of columns
    int cols() {
        return _cols;
    }

    /// Required by Spectra.
    /// Computes the product Cx = y, i.e. @f[ (B^-1 A)v = y@$]. But B = LU, so Ax = LUy.
    /// Let Ax = v, then LUy = v. Then Lw = v and finally Uy = w to get y;
    /// \param[in] x_in
    /// \param[out] y_out
    void perform_op(NT const * x_in, NT* y_out) {

        // Declaring the vectors like this, we don't copy the values of x_in to v
        // and next of y to y_out
        Eigen::Map<VT> const x(const_cast<double*>(x_in), _rows);
        VT const v = *A * x;

        Eigen::Map<VT> y(y_out, _rows);
        y = Blu.solve(v);
    }

    /// Required by arpack.
    /// Computes the product Cx = y, i.e. @f[ (B^-1 A)v = y@$]. But B = LU, so Ax = LUy.
    /// Let Ax = v, then LUy = v. Then Lw = v and finally Uy = w to get y;
    /// \param[in] x_in
    /// \param[out] y_out
    void MultMv(NT * x_in, NT* y_out) {

        // Declaring the vectors like this, we don't copy the values of x_in to v
        // and next of y to y_out
        Eigen::Map<VT> const x(const_cast<double*>(x_in), _rows);
        VT const v = *A * x;

        Eigen::Map<VT> y(y_out, _rows);
        y = Blu.solve(v);
    }
};
#endif //VOLESTI_DENSEPRODUCTMATRIX_H
