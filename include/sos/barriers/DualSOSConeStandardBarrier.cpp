//
// Created by test Bento Natura on 22/07/2020.
//

#include "DualSOSConeStandardBarrier.h"

Vector DualSOSConeStandardBarrier::gradient(Vector x) {
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

Matrix DualSOSConeStandardBarrier::hessian(Vector x) {
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

//FIXME: Wrong implementation in in_interior method.
bool DualSOSConeStandardBarrier::in_interior(Vector x) {
    Matrix X = Lambda(x);
    return X.determinant() > 0;
}

IPMDouble DualSOSConeStandardBarrier::concordance_parameter(Vector x) {
    return x.rows();
}

Vector DualSOSConeStandardBarrier::initialize_x() {
    //TODO: find centered initialization
    return Vector();
}

Vector DualSOSConeStandardBarrier::initialize_s() {
    //TODO: find centered initialization
    return Vector();
}

Matrix DualSOSConeStandardBarrier::Lambda(Vector x) {
    assert(x.rows() == _max_polynomial_degree + 1);
    Matrix M(_max_polynomial_degree + 1, _max_polynomial_degree + 1);
    for (unsigned i = 0; i < _max_polynomial_degree + 1; ++i) {
        for (unsigned j = 0; j < _max_polynomial_degree + 1; ++j) {
            M(i, j) = x(i + j);
        }
    }
    return M;
}

