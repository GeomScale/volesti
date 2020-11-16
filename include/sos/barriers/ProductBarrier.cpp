// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file


#include "ProductBarrier.h"

//TODO: Figure out how to compute llt decomposition with diagonal blocks. There
// does not seem to be a straightforward way with Eigen.
Eigen::LLT<Matrix> ProductBarrier::llt(Vector x, bool symmetrize) {
    return LHSCB::llt(x, symmetrize);
}

//TODO: figure out how these methods can be even more abstracted.

Vector ProductBarrier::llt_L_solve(Vector x, Vector rhs) {
    update_segments();
    Vector product_llt_solve(_num_variables);
#ifdef PARALLELIZE_BARRIERS
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        std::pair<int, int> & seg = _segments[i];
        LHSCB *barrier = _barriers[i];
        Vector x_seg = x.segment(seg.first, seg.second - seg.first);
        Vector rhs_seg = rhs.segment(seg.first, seg.second - seg.first);
        Vector lls_solve_seg = barrier->llt_L_solve(x_seg, rhs_seg);
        product_llt_solve.segment(seg.first , seg.second - seg.first) = lls_solve_seg;
    }
    return product_llt_solve;
}

Matrix ProductBarrier::llt_solve(Vector x, const Matrix &rhs) {
    update_segments();
    Matrix product_llt_solve = Matrix::Zero(rhs.rows(), rhs.cols());
#ifdef PARALLELIZE_BARRIERS
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        std::pair<int, int> & seg = _segments[i];
        LHSCB *barrier = _barriers[i];
        Vector x_seg = x.segment(seg.first, seg.second - seg.first);
        Matrix rhs_block = rhs.block(seg.first, 0, seg.second-seg.first, rhs.cols());
        Matrix lls_solve_block = barrier->llt_solve(x_seg, rhs_block);
        product_llt_solve.block(seg.first, 0, seg.second - seg.first, rhs.cols()) = lls_solve_block;
    }
    return product_llt_solve;
}

//Should not be used, as immense storage is needed.
Matrix ProductBarrier::hessian(Vector x) {
   return evaluate(x, &LHSCB::hessian);
}

bool ProductBarrier::in_interior(Vector x) {
    _in_interior_timer.start();
    update_segments();
    std::vector<bool> in_interior_vec(_barriers.size());
#ifdef PARALLELIZE_BARRIERS
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        std::pair<int, int> & seg = _segments[i];
        LHSCB *barrier = _barriers[i];
        Vector x_seg  = x.segment(seg.first, seg.second - seg.first);
        in_interior_vec[i] =  barrier->in_interior(x_seg);
    }
    _in_interior_timer.stop();
    return std::all_of(in_interior_vec.begin(), in_interior_vec.end(), [](bool b){return b;});
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
    return evaluate(&LHSCB::initialize_s);
}

Vector ProductBarrier::initialize_s() {
    return evaluate(&LHSCB::initialize_s);
}

Matrix ProductBarrier::inverse_hessian(Vector x) {
    return evaluate(x, &LHSCB::inverse_hessian);
}

Vector ProductBarrier::gradient(Vector x) {
    return evaluate(x, &LHSCB::gradient);
}

void ProductBarrier::update_segments() {
    _segments.resize(0);
    unsigned idx = 0;
    std::vector<std::pair<int, int> > segments;
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        unsigned num_variables = _num_vars_per_barrier[i];
        _segments.push_back(std::pair<int, int>(idx, idx + num_variables));
        idx += num_variables;
    }
}

Vector ProductBarrier::evaluate(Vector x, VectorFunc func) {
    update_segments();
    Vector product_vector(_num_variables);
#ifdef PARALLELIZE_BARRIERS
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        std::pair<int, int> &seg = _segments[i];
        LHSCB *barrier = _barriers[i];
        Vector x_seg = x.segment(seg.first, seg.second - seg.first);
        Vector vec_seg = (barrier->*func)(x_seg);
        product_vector.segment(seg.first, seg.second - seg.first) = vec_seg;
    }
    return product_vector;
}

Matrix ProductBarrier::evaluate(Vector x, MatrixFunc func) {
    update_segments();
    Matrix product_matrix = Matrix::Zero(_num_variables, _num_variables);
#ifdef PARALLELIZE_BARRIERS
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        std::pair<int, int> &seg = _segments[i];
        LHSCB *barrier = _barriers[i];
        Vector x_seg = x.segment(seg.first, seg.second - seg.first);
        Matrix matrix_block = (barrier->*func)(x_seg);
        product_matrix.block(seg.first, seg.first, seg.second - seg.first, seg.second - seg.first) = matrix_block;
    }
    return product_matrix;
}

Vector ProductBarrier::evaluate(VoidFunc func) {
    update_segments();
    Vector product_vector(_num_variables);
#ifdef PARALLELIZE_BARRIERS
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < _barriers.size(); ++i) {
        std::pair<int, int> &seg = _segments[i];
        LHSCB *barrier = _barriers[i];
        Vector vec_seg = (barrier->*func)();
        product_vector.segment(seg.first, seg.second - seg.first) = vec_seg;
    }
    return product_vector;
}