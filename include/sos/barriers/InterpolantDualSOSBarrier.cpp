//
// Created by test Bento Natura on 22/07/2020.
//

#include "InterpolantDualSOSBarrier.h"


bool InterpolantDualSOSBarrier::update_gradient_hessian_LLT(Vector x, bool check_interior_only) {
    _preintermediate_matrix.noalias() = _g.cwiseProduct(x).asDiagonal() * _P;
    _intermediate_matrix.noalias() = _P.transpose() * _preintermediate_matrix;
    _intermediate_LLT.compute( _intermediate_matrix);

    if(_intermediate_LLT.info() == Eigen::NumericalIssue){
        return false;
    }

    if(check_interior_only){
        return true;
    }

    _V.noalias() = _intermediate_LLT.matrixL().solve(_P.transpose());
    //Experiments showed that using the triangularView slows the program down.
//    _Q.triangularView<Eigen::Lower>() = _V.transpose() * _V;
//    _Q = _Q.selfadjointView<Eigen::Lower>();
    _Q.noalias() = _V.transpose() * _V;

    //TODO: store hessian as self-adjoint

    if (_stored_hessians.empty()) {
        _stored_hessians.resize(1);
    }
    _stored_hessians[0].first = x;
    _stored_hessians[0].second.noalias() = _g_g_transpose.cwiseProduct(_Q.cwiseProduct(_Q));

    if (_stored_gradients.empty()) {
        _stored_gradients.resize(1);
    }
    _stored_gradients[0].first = x;
    _stored_gradients[0].second.noalias() = -_Q.diagonal().cwiseProduct(_g);

    if (_stored_LLT.empty()) {
        _stored_LLT.resize(1);
    }

    _stored_LLT[0].first = x;
    _stored_LLT[0].second.compute(_stored_hessians[0].second.selfadjointView<Eigen::Lower>());

    return true;
}

Vector InterpolantDualSOSBarrier::gradient(Vector x) {
    auto *grad_ptr = find_gradient(x);
    if (grad_ptr) {
        return *grad_ptr;
    }
    update_gradient_hessian_LLT(x);
    return _stored_gradients[0].second;
}


Matrix InterpolantDualSOSBarrier::hessian(Vector x) {
    auto *hess_ptr = find_hessian(x);
    if (hess_ptr) {
        return *hess_ptr;
    }
    update_gradient_hessian_LLT(x);
    return _stored_hessians[0].second;
}

Eigen::LLT<Matrix> InterpolantDualSOSBarrier::llt(Vector x, bool) {
    auto *llt_ptr = find_LLT(x);
    if (llt_ptr) {
        return *llt_ptr;
    }
    update_gradient_hessian_LLT(x);
    return _stored_LLT[0].second;
}

//Should not be inveoked as it is slow.
Matrix InterpolantDualSOSBarrier::inverse_hessian(Vector x) {
    //Not sure if correct method.
    Matrix L_inv = llt(x).matrixL().toDenseMatrix().inverse();
    Matrix inv = L_inv.transpose() * L_inv;
    return inv;
}


//REMARK: Note that in_interior can return true, while llt or hessian could fail. This is due to numerical instabilities in these more complicated operations.
//Currently fixed by disabling line search.
//TODO: calculate gradient and hessian as we are alaready half way there and need to calculate the gradient for centrality anyway.
bool InterpolantDualSOSBarrier::in_interior(Vector x) {
    return update_gradient_hessian_LLT(x, true);
}

IPMDouble InterpolantDualSOSBarrier::concordance_parameter(Vector) {
    return _L;
}

Vector InterpolantDualSOSBarrier::initialize_x() {
    return Vector::Ones(_U);
}

//TODO: Fix dual initialization
Vector InterpolantDualSOSBarrier::initialize_s() {
    return -gradient(initialize_x());
}

