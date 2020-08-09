// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_INTERPOLANTDUALSOSBARRIER_H
#define NONSYMMETRICCONICOPTIMIZATION_INTERPOLANTDUALSOSBARRIER_H

#include "LHSCB.h"
#include "ChebTools/ChebTools.h"

class InterpolantDualSOSBarrier : public LHSCB {

public:
    InterpolantDualSOSBarrier() : LHSCB() {};

    InterpolantDualSOSBarrier(unsigned max_polynomial_degree_, unsigned num_variables_ = 1) :
            InterpolantDualSOSBarrier(max_polynomial_degree_, Vector::Ones(1), num_variables_) {};

    InterpolantDualSOSBarrier(unsigned max_polynomial_degree_, Vector poly_g, unsigned num_variable_symbols_ = 1);

    bool update_gradient_hessian_LLT(Vector x, bool check_interior_only = false);

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    Eigen::LLT<Matrix> llt(Vector x, bool symmetrize = 0) override;

    Matrix inverse_hessian(Vector x) override;

    bool in_interior(Vector x) override;

    //TODO: better solution for implementation concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

    std::vector<InterpolantDouble> &get_basis() {
        return _unisolvent_basis;
    }

    Matrix get_P() {
        return _P;
    }

private:
    unsigned _max_polynomial_degree;
    unsigned _num_variable_symbols;
    std::vector<InterpolantDouble> _unisolvent_basis;
    Matrix _intermediate_matrix;
    Matrix _preintermediate_matrix;
    Eigen::LLT<Matrix> _intermediate_LLT;
    Matrix _Q;
    Matrix _V;
    unsigned _L, _U;

    void construct_univariate(Vector poly_g);
    void construct_bivariate(Vector poly_g);
    void construct_multivariate(Vector poly_g);

    //Weighted polynomials
    Vector _g;
    Matrix _g_g_transpose;
    Matrix _P;
};


#endif //NONSYMMETRICCONICOPTIMIZATION_INTERPOLANTDUALSOSBARRIER_H
