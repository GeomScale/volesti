// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SOS_FULLSPACEBARRIER_H
#define SOS_FULLSPACEBARRIER_H

#include "LHSCB.h"

template<typename IPMDouble>
class FullSpaceBarrier final : public LHSCB<IPMDouble>{

    using LHSCB = LHSCB<IPMDouble>;

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    Matrix inverse_hessian(Vector x) override;

    bool in_interior(Vector x) override;

    Matrix llt_solve(Vector x, const Matrix &rhs) override;
    Vector llt_L_solve(Vector x, Vector rhs) override;


    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

public:
    FullSpaceBarrier(unsigned num_variables_) {
        this->_num_variables = num_variables_;
    }
};

#include "FullSpaceBarrier.hpp"

#endif //SOS_FULLSPACEBARRIER_H
