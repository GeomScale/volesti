// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file
#ifndef SOS_DUALSOSCONESTANDARDBARRIER_H
#define SOS_DUALSOSCONESTANDARDBARRIER_H

#include "LHSCB.h"

//Implementation for Standard Monomial Basis for the Dual SOS cone.
template <typename IPMDouble>
class DualSOSConeStandardBarrier : public LHSCB<IPMDouble> {

    using LHSCB = LHSCB<IPMDouble>;

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

public:
    DualSOSConeStandardBarrier() : LHSCB() {};

    DualSOSConeStandardBarrier(unsigned max_polynomial_degree_) : _max_polynomial_degree(max_polynomial_degree_) {
    };

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    //Matrix inverse_hessian(Vector x) override;
    bool in_interior(Vector x) override;

    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;
    Vector initialize_s() override;

    //Lambda is the linear operator used in Proposition 1.1 in Papp & Yildiz "SOS
    //without semidefinite programming. In particular, for the monomial basis it takes
    //the form listed in section 3.1 of the same paper.
    Matrix Lambda(Vector x);

private:
    unsigned _max_polynomial_degree;
};

#include "DualSOSConeStandardBarrier.hpp"

#endif //SOS_DUALSOSCONESTANDARDBARRIER_H
