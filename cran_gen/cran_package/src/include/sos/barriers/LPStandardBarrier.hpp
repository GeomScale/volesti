// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "LPStandardBarrier.h"

template<typename IPMDouble>
bool LPStandardBarrier<IPMDouble>::in_interior(Vector x) {
    return (x.minCoeff() > 0);
}

template<typename IPMDouble>
Vector<IPMDouble> LPStandardBarrier<IPMDouble>::gradient(Vector x) {
    assert(in_interior(x));
    return -x.array().inverse();
}

template<typename IPMDouble>
Matrix<IPMDouble> LPStandardBarrier<IPMDouble>::hessian(Vector x) {
    return x.array().pow(2).inverse().matrix().asDiagonal();
}

template<typename IPMDouble>
Matrix<IPMDouble> LPStandardBarrier<IPMDouble>::inverse_hessian(Vector x) {
    return hessian(x).inverse();
}

template<typename IPMDouble>
IPMDouble LPStandardBarrier<IPMDouble>::concordance_parameter(Vector x) {
    return x.rows();
}

template<typename IPMDouble>
Vector<IPMDouble> LPStandardBarrier<IPMDouble>::initialize_x() {
    return Vector::Ones(this->_num_variables);
}

template<typename IPMDouble>
Vector<IPMDouble> LPStandardBarrier<IPMDouble>::initialize_s() {
    return Vector::Ones(this->_num_variables);
}
