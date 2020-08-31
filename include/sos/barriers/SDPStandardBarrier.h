// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_SDPSTANDARDBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_SDPSTANDARDBARRIER_H

#include "LHSCB.h"

template<typename IPMDouble>
class SDPStandardBarrier final : public LHSCB<IPMDouble> {
public:
    using LHSCB = LHSCB<IPMDouble>;

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

    SDPStandardBarrier() : LHSCB() {};

    SDPStandardBarrier(unsigned matrix_dimension_) : _matrix_dimension(matrix_dimension_) {
        this->_num_variables = _matrix_dimension * _matrix_dimension;
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

#include "SDPStandardBarrier.tpp"

#endif //NONSYMMETRICCONICOPTIMIZATION_SDPSTANDARDBARRIER_H
