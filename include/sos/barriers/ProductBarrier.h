//
// Created by test Bento Natura on 22/07/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_PRODUCTBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_PRODUCTBARRIER_H

#include "LHSCB.h"

class ProductBarrier : public LHSCB {

public:
    ProductBarrier() : LHSCB() {}

    ProductBarrier(std::vector<LHSCB *> barriers_, std::vector<unsigned> num_variables_) {
        assert(barriers_.size() == num_variables_.size());
        for (unsigned j = 0; j < barriers_.size(); ++j) {
            _barriers.push_back(barriers_[j]);
            _num_vars_per_barrier.push_back(num_variables_[j]);
        }
    }

    void add_barrier(LHSCB *lhscb) {
        _barriers.push_back(lhscb);
        _num_vars_per_barrier.push_back(lhscb->getNumVariables());
        _num_variables += lhscb->getNumVariables();

    }

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

//    Matrix inverse_hessian(Vector x) override;
    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

    Eigen::LLT<Matrix> llt(Vector x, bool symmetrize = 0) override;

    Matrix llt_solve(Vector x, const Matrix &rhs) override;

    Vector llt_L_solve(Vector x, Vector rhs) override;

    Matrix inverse_hessian(Vector x) override;

private:
    std::vector<LHSCB *> _barriers;
    std::vector<unsigned> _num_vars_per_barrier;
};

//corresponds to intersection of cones.
class SumBarrier : public LHSCB {

public:
    SumBarrier(unsigned num_variables_) : LHSCB() {
        _num_variables = num_variables_;
    }

    SumBarrier(std::vector<LHSCB *> barriers_, unsigned num_variables_) {
        _num_variables = num_variables_;
        for (unsigned j = 0; j < barriers_.size(); ++j) {
            assert(barriers_[j]->getNumVariables() == num_variables_);
            _barriers.push_back(barriers_[j]);
        }
    }

    void add_barrier(LHSCB *lhscb) {
        _barriers.push_back(lhscb);
    }

    Vector gradient(Vector x) override {
        Vector grad_vec = Vector::Zero(x.rows());
        for (auto barrier : _barriers) {
            grad_vec += barrier->gradient(x);
        }
        return grad_vec;
    }

    Matrix hessian(Vector x) override {
        Matrix hess_mat = Matrix::Zero(x.rows(), x.rows());
        for (auto barrier : _barriers) {
            hess_mat += barrier->hessian(x);
        }
        return hess_mat;
    }

    bool in_interior(Vector x) override {
        for (auto barrier : _barriers) {
            if (not barrier->in_interior(x)) {
                return false;
            }
        }
        return true;
    }

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override {
        IPMDouble conc_par = 0.;
        for (auto barrier : _barriers) {
            conc_par += barrier->concordance_parameter(x);
        }
        return conc_par;
    }

    //TODO: verify that initialisation works for both primal and dual side in general.
    Vector initialize_x() override {
        assert(!_barriers.empty());
        return _barriers[0]->initialize_x();
    }

    Vector initialize_s() override {
        Vector s_init = Vector::Zero(_num_variables);
        for (auto barrier : _barriers) {
            s_init += barrier->initialize_s();
        }
        return s_init;
    }

private:
    std::vector<LHSCB *> _barriers;
};



#endif //NONSYMMETRICCONICOPTIMIZATION_PRODUCTBARRIER_H
