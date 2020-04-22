//
// Created by panagiotis on 16/11/2019.
//

#ifndef VOLESTI_SPECTRALU_H
#define VOLESTI_SPECTRALU_H


#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <stdexcept>

namespace Spectra {


    template <typename Scalar, int Uplo = Eigen::Lower, int Flags = 0, typename StorageIndex = int>
    class SpectraLU
    {
    private:
        typedef Eigen::Index Index;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Map<const Vector> MapConstVec;
        typedef Eigen::Map<Vector> MapVec;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

        ConstGenericMatrix m_mat;
        const int m_n;
        Eigen::PartialPivLU<Matrix> m_decomp;

    public:
        ///
        /// Constructor to create the matrix operation object.
        ///
        /// \param mat An **Eigen** sparse matrix object, whose type can be
        /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
        /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
        ///
        SpectraLU(ConstGenericMatrix& mat) :
                m_mat(mat), m_n(mat.rows())
        {
            if(mat.rows() != mat.cols())
                throw std::invalid_argument("SparseRegularInverse: matrix must be square");

            m_decomp.compute(mat);
        }

        ///
        /// Return the number of rows of the underlying matrix.
        ///
        Index rows() const { return m_n; }
        ///
        /// Return the number of columns of the underlying matrix.
        ///
        Index cols() const { return m_n; }

        ///
        /// Perform the solving operation \f$y=B^{-1}x\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = inv(B) * x_in
        void solve(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in,  m_n);
            MapVec      y(y_out, m_n);
            y.noalias() = m_decomp.solve(x);
        }

        ///
        /// Perform the matrix-vector multiplication operation \f$y=Bx\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = B * x_in
        void mat_prod(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in,  m_n);
            MapVec      y(y_out, m_n);
            y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
        }
    };


} // namespace Spectra


//
//namespace Spectra {
//
//
/////
///// \ingroup MatOp
/////
///// This class defines the operations related to Cholesky decomposition on a
///// positive definite matrix, \f$B=LL'\f$, where \f$L\f$ is a lower triangular
///// matrix. It is mainly used in the SymGEigsSolver generalized eigen solver
///// in the Cholesky decomposition mode.
/////
//    template <typename Scalar, int Uplo = Eigen::Lower>
//    class SpectraLU
//    {
//    private:
//        typedef Eigen::Index Index;
//        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
//        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
//        typedef Eigen::Map<const Matrix> MapConstMat;
//        typedef Eigen::Map<const Vector> MapConstVec;
//        typedef Eigen::Map<Vector> MapVec;
//        typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;
//
//        const Index m_n;
//        Eigen::PartialPivLU<Matrix> m_decomp;
//        Matrix L;
//        Matrix U;
//        int m_info;  // status of the decomposition
//
//    public:
//        ///
//        /// Constructor to create the matrix operation object.
//        ///
//        /// \param mat An **Eigen** matrix object, whose type can be
//        /// `Eigen::Matrix<Scalar, ...>` (e.g. `Eigen::MatrixXd` and
//        /// `Eigen::MatrixXf`), or its mapped version
//        /// (e.g. `Eigen::Map<Eigen::MatrixXd>`).
//        ///
//        SpectraLU(ConstGenericMatrix& mat) :
//                m_n(mat.rows()), m_info(NOT_COMPUTED)
//        {
//            if(mat.rows() != mat.cols())
//                throw std::invalid_argument("DenseCholesky: matrix must be square");
//
//            m_decomp.compute(mat);
//
//            int m = mat.rows();
//            U = m_decomp.matrixLU();
//            L = Eigen::MatrixXd::Identity(dim, dim);
//
//            double* data = L.data();
//            double* dataU = U.data();
//
//            for (int at_col = 0; at_col < m-1; at_col++) {
//                int col_offset = at_col * m;
//                double *target = data + col_offset + at_col;
//                double *target = dataU + col_offset + at_col;
//
//                for (int at_row = at_col+1; at_row < m; at_row++) {
//                    *(target++) = *(target++);
//                }
//            }
//
//            data = U.data();
//
//            for (int at_col = 0; at_col < m-1; at_col++) {
//                int col_offset = at_col * m;
//                double *target = data + col_offset + at_col;
//
//                for (int at_row = at_col+1; at_row < m; at_row++) {
//                    *(target++) = 0;
//                }
//            }
//
//            m_info = (m_decomp.info() == Eigen::Success) ?
//                     SUCCESSFUL :
//                     NUMERICAL_ISSUE;
//        }
//
//        ///
//        /// Returns the number of rows of the underlying matrix.
//        ///
//        Index rows() const { return m_n; }
//        ///
//        /// Returns the number of columns of the underlying matrix.
//        ///
//        Index cols() const { return m_n; }
//
//        ///
//        /// Returns the status of the computation.
//        /// The full list of enumeration values can be found in \ref Enumerations.
//        ///
//        int info() const { return m_info; }
//
//        ///
//        /// Performs the lower triangular solving operation \f$y=L^{-1}x\f$.
//        ///
//        /// \param x_in  Pointer to the \f$x\f$ vector.
//        /// \param y_out Pointer to the \f$y\f$ vector.
//        ///
//        // y_out = inv(L) * x_in
//        void lower_triangular_solve(const Scalar* x_in, Scalar* y_out) const
//        {
//            MapConstVec x(x_in,  m_n);
//            MapVec      y(y_out, m_n);
//            y.noalias() = m_decomp.matrixL().solve(x);
//        }
//
//        ///
//        /// Performs the upper triangular solving operation \f$y=(L')^{-1}x\f$.
//        ///
//        /// \param x_in  Pointer to the \f$x\f$ vector.
//        /// \param y_out Pointer to the \f$y\f$ vector.
//        ///
//        // y_out = inv(L') * x_in
//        void upper_triangular_solve(const Scalar* x_in, Scalar* y_out) const
//        {
//            MapConstVec x(x_in,  m_n);
//            MapVec      y(y_out, m_n);
//            y.noalias() = m_decomp.matrixU().solve(x);
//        }
//    };



//} // namespace Spectra


#endif //VOLESTI_SPECTRALU_H
