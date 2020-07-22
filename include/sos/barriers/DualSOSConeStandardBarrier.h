//
// Created by test Bento Natura on 22/07/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_DUALSOSCONESTANDARDBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_DUALSOSCONESTANDARDBARRIER_H

#include "LHSCB.h"

//Implementation for Standard Monomial Basis
class DualSOSConeStandardBarrier : public LHSCB {

public:
    DualSOSConeStandardBarrier() : LHSCB() {};

    DualSOSConeStandardBarrier(unsigned max_polynomial_degree_) : _max_polynomial_degree(max_polynomial_degree_) {
    };

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

//    Matrix inverse_hessian(Vector x) override;
    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

    Matrix Lambda(Vector x);

private:
    unsigned _max_polynomial_degree;
};



#endif //NONSYMMETRICCONICOPTIMIZATION_DUALSOSCONESTANDARDBARRIER_H
