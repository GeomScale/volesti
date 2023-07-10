// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "DualSOSConeStandardBarrier.h"


template <typename IPMDouble>
Vector<IPMDouble> DualSOSConeStandardBarrier<IPMDouble>::gradient(Vector x) {
    assert(in_interior(x));
    Matrix X = Lambda(x);
    Matrix Z = X.inverse();
    Vector g(x.rows());
    for (int i = 0; i < g.rows(); ++i) {
        Matrix E_i = Matrix::Zero(X.rows(), X.cols());
        for (int j = 0; j <= i; ++j) {
            E_i(j, i - j) = 1;
        }
        g(i) = -Z.cwiseProduct(E_i).sum();
    }
    return g;
}

template <typename IPMDouble>
Matrix<IPMDouble> DualSOSConeStandardBarrier<IPMDouble>::hessian(Vector x) {
    assert(in_interior(x));
    Matrix X = Lambda(x);
    Matrix Z = X.inverse();
    Matrix H = Matrix::Zero(Z.rows(), Z.cols());
    for (int u = 0; u < H.rows(); ++u) {
        for (int v = 0; v < H.cols(); ++v) {
            IPMDouble H_uv = 0;
            for (int a = 0; a <= u; ++a) {
                for (int k = 0; k <= v; ++k) {
                    H_uv += Z(a, u - a) + Z(k, v - k);
                }
            }
        }
    }
    return H;
}

template <typename IPMDouble>
bool DualSOSConeStandardBarrier<IPMDouble>::in_interior(Vector x) {
    Matrix X = Lambda(x);
    CustomLLT<Matrix, Eigen::Lower> llt_check;
    llt_check.compute(X);
    return llt_check.info() != Eigen::NumericalIssue;
}

template <typename IPMDouble>
IPMDouble DualSOSConeStandardBarrier<IPMDouble>::concordance_parameter(Vector x) {
    return x.rows();
}

template <typename IPMDouble>
Vector<IPMDouble> DualSOSConeStandardBarrier<IPMDouble>::initialize_x() {
    //TODO: find centered initialization
    return Vector();
}

template <typename IPMDouble>
Vector<IPMDouble> DualSOSConeStandardBarrier<IPMDouble>::initialize_s() {
    //TODO: find centered initialization
    return Vector();
}

template <typename IPMDouble>
Matrix<IPMDouble> DualSOSConeStandardBarrier<IPMDouble>::Lambda(Vector x) {
    assert(x.rows() == _max_polynomial_degree + 1);
    Matrix M(_max_polynomial_degree + 1, _max_polynomial_degree + 1);
    for (unsigned i = 0; i < _max_polynomial_degree + 1; ++i) {
        for (unsigned j = 0; j < _max_polynomial_degree + 1; ++j) {
            M(i, j) = x(i + j);
        }
    }
    return M;
}

