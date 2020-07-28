// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "FullSpaceBarrier.h"

//This barrier function exists for convenience.
//An instance should be reformulated instead of using the barrier.

Vector FullSpaceBarrier::gradient(Vector) {
    return Vector::Zero(_num_variables);
}

Matrix FullSpaceBarrier::hessian(Vector) {
    return Matrix::Zero(_num_variables, _num_variables);
}

Matrix FullSpaceBarrier::inverse_hessian(Vector) {
    return Matrix::Zero(_num_variables, _num_variables);
}

bool FullSpaceBarrier::in_interior(Vector) {
    return true;
}

Matrix FullSpaceBarrier::llt_solve(Vector x, const Matrix &rhs) {
    if(rhs.norm() > 1e-10){
        spdlog::warn("Exit because RHS of matrix is");
        std::cout << rhs << std::endl;
        assert(false);
    }
    return Matrix::Zero(rhs.rows(),rhs.cols());
}

Vector FullSpaceBarrier::llt_L_solve(Vector x, Vector rhs) {
    if(rhs.norm() > 1e-10){
        spdlog::warn("Exit because RHS of matrix is");
        std::cout << rhs << std::endl;
        assert(false);
    }
    return Vector::Zero(x.rows());
}


IPMDouble FullSpaceBarrier::concordance_parameter(Vector) {
    return 0;
}

Vector FullSpaceBarrier::initialize_x() {
    return Vector::Zero(_num_variables);
}

Vector FullSpaceBarrier::initialize_s() {
    return Vector::Zero(_num_variables);
}
