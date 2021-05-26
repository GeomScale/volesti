//
// Created by test Bento Natura on 30/07/2020.
//

#ifndef SOS_SUMBARRIER_H
#define SOS_SUMBARRIER_H

#include "LHSCB.h"
#include "InterpolantDualSOSBarrier.h"

//corresponds to intersection of cones.
template<typename IPMDouble>
class SumBarrier : public LHSCB<IPMDouble> {

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

public:
    SumBarrier(unsigned num_variables_);

    SumBarrier(std::vector<LHSCB<IPMDouble> *> barriers_, unsigned num_variables_);

    template<typename T>
    SumBarrier<T> * cast(){
        SumBarrier<T> * sb = new SumBarrier<T>(this->_num_variables);
        for(LHSCB<IPMDouble> * barrier : _barriers){
            //TODO: High priority. Figure out how to undo the cast to the InterpolantBarrier.
            InterpolantDualSOSBarrier<T> * new_barr = static_cast<InterpolantDualSOSBarrier<IPMDouble>*>(barrier)->template cast<T>();
//            LHSCB<T> * new_barr = (*barrier).template cast<T>();
            sb->add_barrier(new_barr);
        }
        return sb;
    };

    void add_barrier(LHSCB<IPMDouble> *lhscb);

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    //TODO: verify that initialisation works for both primal and dual side in general.
    Vector initialize_x() override;

    Vector initialize_s() override;

    std::vector<LHSCB<IPMDouble> *> & get_barriers(){
        return _barriers;
    }
private:
    std::vector<LHSCB<IPMDouble> *> _barriers;
};

#include "SumBarrier.hpp"
#endif //SOS_SUMBARRIER_H
