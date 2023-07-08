// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SOS_LPSTANDARDBARRIER_H
#define SOS_LPSTANDARDBARRIER_H

#include "LHSCB.h"

template <typename IPMDouble>
class LPStandardBarrier final : public LHSCB<IPMDouble> {
    using LHSCB = LHSCB<IPMDouble>;

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

public:
    LPStandardBarrier() : LHSCB() {};

    LPStandardBarrier(unsigned num_variables_) {
        this->_num_variables = num_variables_;
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

#include "LPStandardBarrier.hpp"

#endif //SOS_LPSTANDARDBARRIER_H
