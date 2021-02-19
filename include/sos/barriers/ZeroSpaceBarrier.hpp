// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "ZeroSpaceBarrier.h"

template <typename IPMDouble>
Vector<IPMDouble> ZeroSpaceBarrier<IPMDouble>::gradient(Vector x) {
    //should not be used
    assert(false);
    return -std::numeric_limits<IPMDouble>::infinity() * Vector::Ones(x.cols());
}

template <typename IPMDouble>
Matrix<IPMDouble> ZeroSpaceBarrier<IPMDouble>::hessian(Vector x) {
    //should not be used
    assert(false);
    return std::numeric_limits<IPMDouble>::infinity() * Matrix::Identity(x.cols(), x.cols());
}

template <typename IPMDouble>
Matrix<IPMDouble> ZeroSpaceBarrier<IPMDouble>::inverse_hessian(Vector) {
    return Matrix::Zero(this->_num_variables, this->_num_variables);
}

template <typename IPMDouble>
bool ZeroSpaceBarrier<IPMDouble>::in_interior(Vector) {
    return true;
}

template <typename IPMDouble>
IPMDouble ZeroSpaceBarrier<IPMDouble>::concordance_parameter(Vector) {
    return 0;
}

template <typename IPMDouble>
Vector<IPMDouble> ZeroSpaceBarrier<IPMDouble>::initialize_x() {
    return Vector::Zero(this->_num_variables);
}

template <typename IPMDouble>
Vector<IPMDouble> ZeroSpaceBarrier<IPMDouble>::initialize_s() {
    return Vector::Zero(this->_num_variables);
}
