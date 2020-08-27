//
// Created by test Bento Natura on 30/07/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_SUMBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_SUMBARRIER_H

#include "LHSCB.h"

//corresponds to intersection of cones.
class SumBarrier : public LHSCB {

public:
    SumBarrier(unsigned num_variables_);

    SumBarrier(std::vector<LHSCB *> barriers_, unsigned num_variables_);

    void add_barrier(LHSCB *lhscb);

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    //TODO: verify that initialisation works for both primal and dual side in general.
    Vector initialize_x() override;

    Vector initialize_s() override;

    std::vector<LHSCB *> & get_barriers(){
        return _barriers;
    }
private:
    std::vector<LHSCB *> _barriers;
};
#endif //NONSYMMETRICCONICOPTIMIZATION_SUMBARRIER_H
