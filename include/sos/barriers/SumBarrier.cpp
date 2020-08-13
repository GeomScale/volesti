//
// Created by test Bento Natura on 30/07/2020.
//

#include "SumBarrier.h"

SumBarrier::SumBarrier(unsigned num_variables_) : LHSCB() {
    _num_variables = num_variables_;
}

SumBarrier::SumBarrier(std::vector<LHSCB *> barriers_, unsigned num_variables_) {
    _num_variables = num_variables_;
    for (unsigned j = 0; j < barriers_.size(); ++j) {
        assert(barriers_[j]->getNumVariables() == num_variables_);
        _barriers.push_back(barriers_[j]);
    }
}

void SumBarrier::add_barrier(LHSCB *lhscb) {
    _barriers.push_back(lhscb);
}

Vector SumBarrier::gradient(Vector x) {
    Vector grad_vec = Vector::Zero(x.rows());
    std::vector<Vector> barrier_vec(_barriers.size());
#pragma omp parallel for
    for (int i = 0; i < _barriers.size(); i++){
        barrier_vec[i] = _barriers[i]->gradient(x);
    }
    for(int i = 0; i < _barriers.size(); i++){
        grad_vec += barrier_vec[i];
    }
    return grad_vec;
}

Matrix SumBarrier::hessian(Vector x) {
    Matrix hess_mat = Matrix::Zero(x.rows(), x.rows());
    std::vector<Matrix> barrier_vec(_barriers.size());
#pragma omp parallel for
    for (int i = 0; i < _barriers.size(); i++){
        barrier_vec[i] = _barriers[i]->hessian(x);
    }
    for(int i = 0; i < _barriers.size(); i++){
        hess_mat += barrier_vec[i];
    }

    return hess_mat;
}

bool SumBarrier::in_interior(Vector x) {
    std::vector<bool> barrier_vec(_barriers.size());
#pragma omp parallel for
    for (int i = 0; i < _barriers.size(); i++){
        barrier_vec[i] = _barriers[i]->in_interior(x);
    }
    return std::all_of(barrier_vec.begin(), barrier_vec.end(), [](bool b){return b;});
}

//TODO: better solution for concordance parameter;
IPMDouble SumBarrier::concordance_parameter(Vector x) {
    IPMDouble conc_par = 0.;
    for (
        auto barrier
            : _barriers) {
        conc_par += barrier->
                concordance_parameter(x);
    }
    return conc_par;
}

//TODO: verify that initialisation works for both primal and dual side in general.
Vector SumBarrier::initialize_x() {
    assert(!_barriers.empty());
    return _barriers[0]->initialize_x();

}

Vector SumBarrier::initialize_s() {
    Vector s_init = Vector::Zero(_num_variables);
    for (
        auto barrier
            : _barriers) {
        s_init += barrier->initialize_s();
    }
    return s_init;
}

