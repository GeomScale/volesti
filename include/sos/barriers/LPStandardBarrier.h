//
// Created by test Bento Natura on 22/07/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_LPSTANDARDBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_LPSTANDARDBARRIER_H

#include "LHSCB.h"

class LPStandardBarrier final : public LHSCB {
public:
    LPStandardBarrier() : LHSCB() {};

    LPStandardBarrier(unsigned num_variables_) {
        _num_variables = num_variables_;
    };

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    Matrix inverse_hessian(Vector x) override;

    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

private:

};

#endif //NONSYMMETRICCONICOPTIMIZATION_LPSTANDARDBARRIER_H