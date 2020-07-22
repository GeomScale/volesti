// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include<iostream>
#include "LHSCB.h"

Matrix LHSCB::inverse_hessian(Vector x) {
    Eigen::LLT<Matrix> LLT = llt(x);
    Matrix L_inv = LLT.matrixL().toDenseMatrix().inverse();
    return L_inv.transpose() * L_inv;
}

unsigned int LHSCB::getNumVariables() const {
    return _num_variables;
}

//TODO: use short queue instead.
Vector *LHSCB::find_gradient(Vector x) {
    for (int i = _stored_gradients.size() - 1; i >= 0; i--) {
        if (x == _stored_gradients[i].first) {
            return &_stored_gradients[i].second;
        }
    }
    return nullptr;
}

Matrix *LHSCB::find_hessian(Vector x) {
    for (int i = _stored_hessians.size() - 1; i >= 0; i--) {
        if (x == _stored_hessians[i].first) {
            return &_stored_hessians[i].second;
        }
    }
    return nullptr;
}

Eigen::LLT<Matrix> *LHSCB::find_LLT(Vector x) {
    for (int i = _stored_LLT.size() - 1; i >= 0; i--) {
        if (x == _stored_LLT[i].first) {
            return &_stored_LLT[i].second;
        }
    }
    return nullptr;
}

