//
// Created by test Bento Natura on 22/07/2020.
//

#include "FullSpaceBarrier.h"

//Should not be used. Reformulate instance instead.

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

IPMDouble FullSpaceBarrier::concordance_parameter(Vector) {
    return 0;
}

Vector FullSpaceBarrier::initialize_x() {
    return Vector::Zero(_num_variables);
}

Vector FullSpaceBarrier::initialize_s() {
    return Vector::Zero(_num_variables);
}
