//
// Created by test Bento Natura on 22/07/2020.
//

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

    const Matrix &inverse_hessian = x.array().pow(2).matrix().asDiagonal();
//    std::cout << "Inverse hessian: " << std::endl << inverse_hessian << std::endl;
    return inverse_hessian;
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
