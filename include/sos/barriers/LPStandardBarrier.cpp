// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "LPStandardBarrier.h"

bool LPStandardBarrier::in_interior(Vector x) {
    return (x.minCoeff() > 0);
}

Vector LPStandardBarrier::gradient(Vector x) {
    assert(in_interior(x));
    return -x.array().inverse();
}

Matrix LPStandardBarrier::hessian(Vector x) {
    return x.array().pow(2).inverse().matrix().asDiagonal();
}

Matrix LPStandardBarrier::inverse_hessian(Vector x) {
    return hessian(x).inverse();
}

IPMDouble LPStandardBarrier::concordance_parameter(Vector x) {
    return x.rows();
}

Vector LPStandardBarrier::initialize_x() {
    return Vector::Ones(_num_variables);
}

Vector LPStandardBarrier::initialize_s() {
    return Vector::Ones(_num_variables);
}
