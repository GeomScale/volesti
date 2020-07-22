//
// Created by test Bento Natura on 22/07/2020.
//

#include "SDPStandardBarrier.h"

bool SDPStandardBarrier::in_interior(Vector x) {
    Matrix X = toMatrix(x);
    auto LLT = X.llt();
    if (LLT.info() != Eigen::NumericalIssue) {
        return true;
    }
    return false;
}

Vector SDPStandardBarrier::gradient(Vector x) {
    assert(in_interior(x));
    Matrix X = toMatrix(x);
    Matrix X_Inv = X.inverse();
    Vector x_inv = toVector(X_Inv);
    return -x_inv;
}

Matrix SDPStandardBarrier::hessian(Vector x) {
    assert(in_interior(x));
    Matrix X = toMatrix(x);
    Matrix X_Inv = X.inverse();

    Matrix H(_matrix_dimension * _matrix_dimension, _matrix_dimension * _matrix_dimension);
    for (unsigned i = 0; i < _matrix_dimension; i++) {
        for (unsigned j = 0; j < _matrix_dimension; ++j) {
            Matrix x_ij = X_Inv.col(j) * X_Inv.row(i);
            Eigen::Map<Matrix> x_ij_row(x_ij.data(), _matrix_dimension * _matrix_dimension, 1);
            H.col(i * _matrix_dimension + j) = x_ij_row;
        }
    }
    return H;
}

//TODO: figure out what correct concordance for SDP is.
IPMDouble SDPStandardBarrier::concordance_parameter(Vector) {
    return _matrix_dimension;
}

Vector SDPStandardBarrier::toVector(Matrix X) {
    assert(X.rows() == _matrix_dimension and X.cols() == _matrix_dimension);
    return MatrixToVector(X);
}

Matrix SDPStandardBarrier::toMatrix(Vector x) {
    assert(x.rows() == _matrix_dimension * _matrix_dimension);
    return VectorToSquareMatrix(x, _matrix_dimension);
}

Vector SDPStandardBarrier::initialize_x() {
    Matrix Id = Matrix::Identity(_matrix_dimension, _matrix_dimension);
    Eigen::Map<Matrix> Id_stacked(Id.data(), _matrix_dimension * _matrix_dimension, 1);
    return Id_stacked;
}

Vector SDPStandardBarrier::initialize_s() {
    return initialize_x();
}

//TODO: Nicer implementation of all the product operations.



