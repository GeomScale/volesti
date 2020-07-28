// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_FULLSPACEBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_FULLSPACEBARRIER_H

#include "LHSCB.h"

class FullSpaceBarrier final : public LHSCB {

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
        _num_variables = num_variables_;
    }
};

#endif //NONSYMMETRICCONICOPTIMIZATION_FULLSPACEBARRIER_H
