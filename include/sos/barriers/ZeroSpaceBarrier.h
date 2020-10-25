// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SOS_ZEROSPACEBARRIER_H
#define SOS_ZEROSPACEBARRIER_H

#include "LHSCB.h"

//This class models the 0-cone. Should not be used, so reformulate your
//instance.
template <typename IPMDouble>
class ZeroSpaceBarrier final : public LHSCB<IPMDouble> {

    using LHSCB = LHSCB<IPMDouble>;

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    Matrix inverse_hessian(Vector x) override;

    bool in_interior(Vector x) override;

    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

public:
    ZeroSpaceBarrier(unsigned num_variables_) {
        this->_num_variables = num_variables_;
    }
};

#include "ZeroSpaceBarrier.hpp"

#endif //SOS_ZEROSPACEBARRIER_H
