//
// Created by test Bento Natura on 22/07/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_SDPSTANDARDBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_SDPSTANDARDBARRIER_H

#include "LHSCB.h"

class SDPStandardBarrier final : public LHSCB {
public:

    SDPStandardBarrier() : LHSCB() {};

    SDPStandardBarrier(unsigned matrix_dimension_) : _matrix_dimension(matrix_dimension_) {
        _num_variables = _matrix_dimension * _matrix_dimension;
    };

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

//    Matrix inverse_hessian(Vector x) override;
    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;


    Matrix toMatrix(Vector x);

    Vector toVector(Matrix X);

private:
    unsigned _matrix_dimension;
};



#endif //NONSYMMETRICCONICOPTIMIZATION_SDPSTANDARDBARRIER_H
