// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SOS_PRODUCTBARRIER_H
#define SOS_PRODUCTBARRIER_H

#include "LHSCB.h"
#include "SumBarrier.h"
#include "InterpolantDualSOSBarrier.h"


template<typename IPMDouble>
class ProductBarrier : public LHSCB<IPMDouble> {

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

    typedef Vector (LHSCB<IPMDouble>::*VectorFunc)(Vector);
    typedef Matrix (LHSCB<IPMDouble>::*MatrixFunc)(Vector);
    typedef Vector (LHSCB<IPMDouble>::*VoidFunc)();

public:
    ProductBarrier(unsigned num_threads = 1) : LHSCB<IPMDouble>() {
        _num_threads = num_threads;
    }

    ProductBarrier(std::vector<LHSCB<IPMDouble> *> barriers_, std::vector<unsigned> num_variables_, unsigned num_threads = 1) {
        assert(barriers_.size() == num_variables_.size());
        _num_threads = num_threads;
        for (unsigned j = 0; j < barriers_.size(); ++j) {
            _barriers.push_back(barriers_[j]);
            _num_vars_per_barrier.push_back(num_variables_[j]);
        }
    }

    template<typename T>
    ProductBarrier<T> *  cast(){
        ProductBarrier<T> *  pb = new ProductBarrier<T>();
        pb->_num_threads = _num_threads;
        pb->_segments = _segments;
        for(LHSCB<IPMDouble> * barrier : _barriers){
            //TODO: High priority. Figure out how to undo the cast to the SumBarrier.
            LHSCB<T> * new_barr;
            if(dynamic_cast<InterpolantDualSOSBarrier<IPMDouble>*>(barrier)){
                new_barr = static_cast<InterpolantDualSOSBarrier<IPMDouble>*>(barrier)->template cast<T>();
            } else {
               if(dynamic_cast<SumBarrier<IPMDouble>*>(barrier)){
                   new_barr = static_cast<SumBarrier<IPMDouble>*>(barrier)->template cast<T>();
               } else {
                   exit(1);
               }
            }
            pb->add_barrier(new_barr);
        }
        return pb;
    };

    void add_barrier(LHSCB<IPMDouble> *lhscb) {
        _barriers.push_back(lhscb);
        _num_vars_per_barrier.push_back(lhscb->getNumVariables());
        this->_num_variables += lhscb->getNumVariables();
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
    std::vector<std::pair<int ,int> > & get_segments(){
        return _segments;
    };

    Vector evaluate(Vector x, VectorFunc func);
    Matrix evaluate(Vector x, MatrixFunc func);
    Vector evaluate(VoidFunc func);

    std::vector<LHSCB<IPMDouble> *> & get_barriers(){
        return _barriers;
    }

    std::vector<LHSCB<IPMDouble> *> _barriers;
    std::vector<std::pair<int ,int> > _segments;
    std::vector<unsigned> _num_vars_per_barrier;
    unsigned _num_threads;

private:

};

#include "ProductBarrier.hpp"

#endif //SOS_PRODUCTBARRIER_H
