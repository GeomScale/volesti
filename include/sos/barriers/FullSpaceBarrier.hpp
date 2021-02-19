// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "FullSpaceBarrier.h"

//This barrier function exists for convenience.
//An instance should be reformulated instead of using the barrier.

template<typename IPMDouble>
Vector<IPMDouble> FullSpaceBarrier<IPMDouble>::gradient(Vector) {
    return Vector::Zero(this->_num_variables);
}

template<typename IPMDouble>
Matrix<IPMDouble> FullSpaceBarrier<IPMDouble>::hessian(Vector) {
    return Matrix::Zero(this->_num_variables, this->_num_variables);
}

template<typename IPMDouble>
Matrix<IPMDouble> FullSpaceBarrier<IPMDouble>::inverse_hessian(Vector) {
    return Matrix::Zero(this->_num_variables, this->_num_variables);
}

template<typename IPMDouble>
bool FullSpaceBarrier<IPMDouble>::in_interior(Vector) {
    return true;
}

template<typename IPMDouble>
Matrix<IPMDouble> FullSpaceBarrier<IPMDouble>::llt_solve(Vector x, const Matrix &rhs) {
    if(rhs.norm() > 1e-10){
        spdlog::warn("Exit because RHS of matrix is");
        std::cout << rhs << std::endl;
        assert(false);
    }
    return Matrix::Zero(rhs.rows(),rhs.cols());
}

template<typename IPMDouble>
Vector<IPMDouble> FullSpaceBarrier<IPMDouble>::llt_L_solve(Vector x, Vector rhs) {
    if(rhs.norm() > 1e-10){
        spdlog::warn("Exit because RHS of matrix is");
        std::cout << rhs << std::endl;
        assert(false);
    }
    return Vector::Zero(x.rows());
}


template<typename IPMDouble>
IPMDouble FullSpaceBarrier<IPMDouble>::concordance_parameter(Vector) {
    return 0;
}

template<typename IPMDouble>
Vector<IPMDouble> FullSpaceBarrier<IPMDouble>::initialize_x() {
    return Vector::Zero(this->_num_variables);
}

template<typename IPMDouble>
Vector<IPMDouble> FullSpaceBarrier<IPMDouble>::initialize_s() {
    return Vector::Zero(this->_num_variables);
}
