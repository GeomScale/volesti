//
// Created by test Bento Natura on 22/07/2020.
//

#include "ZeroSpaceBarrier.h"


Vector ZeroSpaceBarrier::gradient(Vector x) {
    //should not be used
    assert(false);
    return -std::numeric_limits<IPMDouble>::infinity() * Vector::Ones(x.cols());
}

Matrix ZeroSpaceBarrier::hessian(Vector x) {
    //should not be used
    assert(false);
    return std::numeric_limits<IPMDouble>::infinity() * Matrix::Identity(x.cols(), x.cols());
}

Matrix ZeroSpaceBarrier::inverse_hessian(Vector) {
    return Matrix::Zero(_num_variables, _num_variables);
}

bool ZeroSpaceBarrier::in_interior(Vector) {
    return true;
}

IPMDouble ZeroSpaceBarrier::concordance_parameter(Vector) {
    return 0;
}

Vector ZeroSpaceBarrier::initialize_x() {
    return Vector::Zero(_num_variables);
}

Vector ZeroSpaceBarrier::initialize_s() {
    return Vector::Zero(_num_variables);
}
