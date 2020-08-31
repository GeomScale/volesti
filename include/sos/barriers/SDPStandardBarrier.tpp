// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "SDPStandardBarrier.h"

template<typename IPMDouble>
bool SDPStandardBarrier<IPMDouble>::in_interior(Vector x) {
    Matrix X = toMatrix(x);
    auto LLT = X.llt();
    if (LLT.info() != Eigen::NumericalIssue) {
        return true;
    }
    std::cout << "Try symmetrizing" << std::endl;
    LLT = ((X + X.transpose())/2).llt();
    if (LLT.info() != Eigen::NumericalIssue) {
        return true;
    }
    return false;
}

template<typename IPMDouble>
Vector<IPMDouble> SDPStandardBarrier<IPMDouble>::gradient(Vector x) {
    assert(in_interior(x));
    Matrix X = toMatrix(x);
    Matrix X_Inv = X.inverse();
    Vector x_inv = toVector(X_Inv);
    return -x_inv;
}

template<typename IPMDouble>
Matrix<IPMDouble> SDPStandardBarrier<IPMDouble>::hessian(Vector x) {
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
template<typename IPMDouble>
IPMDouble SDPStandardBarrier<IPMDouble>::concordance_parameter(Vector) {
    return _matrix_dimension;
}

template<typename IPMDouble>
Vector<IPMDouble> SDPStandardBarrier<IPMDouble>::toVector(Matrix X) {
    assert(X.rows() == _matrix_dimension and X.cols() == _matrix_dimension);
    return StackMatrixToVector(X);
}

template<typename IPMDouble>
Matrix<IPMDouble> SDPStandardBarrier<IPMDouble>::toMatrix(Vector x) {
    assert(x.rows() == _matrix_dimension * _matrix_dimension);
    return UnstackVectorToMatrix(x, _matrix_dimension);
}

template<typename IPMDouble>
Vector<IPMDouble> SDPStandardBarrier<IPMDouble>::initialize_x() {
    Matrix Id = Matrix::Identity(_matrix_dimension, _matrix_dimension);
    Eigen::Map<Matrix> Id_stacked(Id.data(), _matrix_dimension * _matrix_dimension, 1);
    return Id_stacked;
}

template<typename IPMDouble>
Vector<IPMDouble> SDPStandardBarrier<IPMDouble>::initialize_s() {
    return initialize_x();
}




