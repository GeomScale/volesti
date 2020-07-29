// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "ProductBarrier.h"

Eigen::LLT<Matrix> ProductBarrier::llt(Vector x, bool symmetrize) {
    return LHSCB::llt(x, symmetrize);
}

//TODO: code resembles code for gradient and other methods. Find abstraction.
Vector ProductBarrier::llt_L_solve(Vector x, Vector rhs) {
    unsigned idx = 0;
    Vector product_llt_solve(_num_variables);
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector x_seg = x.segment(idx, num_variables);
        Vector rhs_seg = rhs.segment(idx, num_variables);
        Vector lls_solve_seg = barrier->llt_L_solve(x_seg, rhs_seg);
        product_llt_solve.segment(idx, num_variables) = lls_solve_seg;
        idx += num_variables;
    }
    return product_llt_solve;
}

//TODO: code resembles code for gradient and other methods. Find abstraction.
Matrix ProductBarrier::llt_solve(Vector x, const Matrix &rhs) {
    unsigned idx = 0;
    Matrix product_llt_solve = Matrix::Zero(rhs.rows(), rhs.cols());
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector x_seg = x.segment(idx, num_variables);
        Matrix rhs_block = rhs.block(idx, 0, num_variables, rhs.cols());
        Matrix lls_solve_block = barrier->llt_solve(x_seg, rhs_block);
        product_llt_solve.block(idx, 0, num_variables, rhs.cols()) = lls_solve_block;
        idx += num_variables;
    }
    return product_llt_solve;
}

//Should not be used, as immense storage is needed.
Matrix ProductBarrier::hessian(Vector x) {
    unsigned idx = 0;
    Matrix product_hessian = Matrix::Zero(_num_variables, _num_variables);
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector v_segment = x.segment(idx, num_variables);
        Matrix hessian_block = barrier->hessian(v_segment);
        product_hessian.block(idx, idx, num_variables, num_variables) = hessian_block;
        idx += num_variables;
    }
    return product_hessian;
}

bool ProductBarrier::in_interior(Vector x) {
    _in_interior_timer.start();
    unsigned idx = 0;
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector v_segment = x.segment(idx, num_variables);
        if (not barrier->in_interior(v_segment)) {
            _in_interior_timer.stop();
            return false;
        }
        idx += num_variables;
    }
    _in_interior_timer.stop();
    return true;
}

IPMDouble ProductBarrier::concordance_parameter(Vector x) {
    unsigned idx = 0;
    IPMDouble concordance_par = 0;
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector v_segment = x.segment(idx, num_variables);
        concordance_par += barrier->concordance_parameter(v_segment);
        idx += num_variables;
    }
    return concordance_par;
}

Vector ProductBarrier::initialize_x() {
    unsigned idx = 0;
    Vector product_init = Vector(_num_variables);
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector init_segment = barrier->initialize_x();
        product_init.segment(idx, num_variables) = init_segment;
        idx += num_variables;
    }
    return product_init;
}

Vector ProductBarrier::initialize_s() {
    unsigned idx = 0;
    Vector product_init = Vector(_num_variables);
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector init_segment = barrier->initialize_s();
        product_init.segment(idx, num_variables) = init_segment;
        idx += num_variables;
    }
    return product_init;
}

Matrix ProductBarrier::inverse_hessian(Vector x) {
    unsigned idx = 0;
    Matrix product_inverse_hessian = Matrix::Zero(_num_variables, _num_variables);
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector v_segment = x.segment(idx, num_variables);
        Matrix inverse_hessian_block = barrier->inverse_hessian(v_segment);
        product_inverse_hessian.block(idx, idx, num_variables, num_variables) = inverse_hessian_block;
        idx += num_variables;
    }
    return product_inverse_hessian;
}

Vector ProductBarrier::gradient(Vector x) {
    unsigned idx = 0;
    Vector product_gradient(_num_variables);
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        LHSCB *barrier = _barriers[i];
        unsigned num_variables = _num_vars_per_barrier[i];
        Vector v_segment = x.segment(idx, num_variables);
        Vector gradient_segment = barrier->gradient(v_segment);
        product_gradient.segment(idx, num_variables) = gradient_segment;
        idx += num_variables;
    }
    return product_gradient;
}