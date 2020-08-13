// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_PRODUCTBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_PRODUCTBARRIER_H

#include "LHSCB.h"

typedef Vector (LHSCB::*VectorFunc)(Vector);
typedef Matrix (LHSCB::*MatrixFunc)(Vector);
typedef Vector (LHSCB::*VoidFunc)();

class ProductBarrier : public LHSCB {

public:
    ProductBarrier(unsigned num_threads = 1) : LHSCB() {
        _num_threads = num_threads;
    }

    ProductBarrier(std::vector<LHSCB *> barriers_, std::vector<unsigned> num_variables_, unsigned num_threads = 1) {
        assert(barriers_.size() == num_variables_.size());
        _num_threads = num_threads;
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

    void update_segments();

    Vector evaluate(Vector x, VectorFunc func);
    Matrix evaluate(Vector x, MatrixFunc func);
    Vector evaluate(VoidFunc func);

private:
    std::vector<LHSCB *> _barriers;
    std::vector<std::pair<int ,int> > _segments;
    std::vector<unsigned> _num_vars_per_barrier;
    unsigned _num_threads;
};

#endif //NONSYMMETRICCONICOPTIMIZATION_PRODUCTBARRIER_H
