// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include<iostream>
#include "LHSCB.h"

template <typename IPMDouble>
Matrix<IPMDouble> LHSCB<IPMDouble>::inverse_hessian(Vector x) {
    Eigen::LLT<Matrix> LLT = llt(x);
    Matrix L_inv = LLT.matrixL().toDenseMatrix().inverse();
    return L_inv.transpose() * L_inv;
}

template <typename IPMDouble>
unsigned int LHSCB<IPMDouble>::getNumVariables() const {
    return _num_variables;
}

//TODO: use short queue instead.
template <typename IPMDouble>
Vector<IPMDouble> *LHSCB<IPMDouble>::find_gradient(Vector x) {
    for (int i = _stored_gradients.size() - 1; i >= 0; i--) {
        if (x == _stored_gradients[i].first) {
            return &_stored_gradients[i].second;
        }
    }
    return nullptr;
}

template <typename IPMDouble>
Matrix<IPMDouble> *LHSCB<IPMDouble>::find_hessian(Vector x) {
    for (int i = _stored_hessians.size() - 1; i >= 0; i--) {
        if (x == _stored_hessians[i].first) {
            return &_stored_hessians[i].second;
        }
    }
    return nullptr;
}

template <typename IPMDouble>
Eigen::LLT<Matrix<IPMDouble> > *LHSCB<IPMDouble>::find_LLT(Vector x) {
    for (int i = _stored_LLT.size() - 1; i >= 0; i--) {
        if (x == _stored_LLT[i].first) {
            return &_stored_LLT[i].second;
        }
    }
    return nullptr;
}

template <typename IPMDouble>
Eigen::LLT<Matrix<IPMDouble> > LHSCB<IPMDouble>::llt(Vector x, bool symmetrize) {

    Eigen::LLT<Matrix> *llt_ptr = nullptr;

    if (not symmetrize) {
        llt_ptr = find_LLT(x);
    }

    if (llt_ptr) {
        return *llt_ptr;
    }

    if (_stored_LLT.empty()) {
        _stored_LLT.resize(1);
    }

    _stored_LLT[0].first = x;
    if (not symmetrize) {
        _stored_LLT[0].second = Eigen::LLT<Matrix>(hessian(x).llt());
    } else {
        Matrix hess_tmp = hessian(x);
        _stored_LLT[0].second = Eigen::LLT<Matrix>(((hess_tmp+ hess_tmp.transpose()) / 2).llt());
    }

    return _stored_LLT[0].second;
}

template <typename IPMDouble>
Matrix<IPMDouble> LHSCB<IPMDouble>::llt_solve(Vector x, const Matrix &rhs) {
    return llt(x).solve(rhs);
}

template <typename IPMDouble>
Vector<IPMDouble> LHSCB<IPMDouble>::llt_L_solve(Vector x, Vector rhs) {
    return llt(x).matrixL().solve(rhs);
}
