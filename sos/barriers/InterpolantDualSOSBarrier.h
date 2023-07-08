// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SOS_INTERPOLANTDUALSOSBARRIER_H
#define SOS_INTERPOLANTDUALSOSBARRIER_H

#include "LHSCB.h"

template <typename IPMDouble>
class InterpolantDualSOSBarrier : public LHSCB<IPMDouble> {

    using LHSCB = LHSCB<IPMDouble>;

    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;

public:
    InterpolantDualSOSBarrier() : LHSCB() {};

    InterpolantDualSOSBarrier(unsigned max_polynomial_degree_, unsigned num_variables_ = 1) :
            InterpolantDualSOSBarrier(max_polynomial_degree_, Vector::Ones(1), num_variables_) {};

    InterpolantDualSOSBarrier(unsigned max_polynomial_degree_, Vector poly_g, unsigned num_variable_symbols_ = 1);

    template<typename T>
    InterpolantDualSOSBarrier<T> * cast(){
        InterpolantDualSOSBarrier<T> * int_bar =  new InterpolantDualSOSBarrier<T>();
        int_bar->_custom_timers.resize(10);
        int_bar->_max_polynomial_degree = _max_polynomial_degree;
        int_bar->_num_variable_symbols = _num_variable_symbols;
        int_bar->_unisolvent_basis = _unisolvent_basis;
        int_bar->_intermediate_matrix = _intermediate_matrix.template cast<T>();
        int_bar->_preintermediate_matrix = _preintermediate_matrix.template cast<T>();
//        int_bar->_intermediate_LLT = _intermediate_LLT.template cast<T>();
        int_bar->_Q = _Q.template cast<T>();
        int_bar->_V = _V.template cast<T>();
        int_bar->_L = _L;
        int_bar->_U = _U;
        int_bar->_g = _g.template cast<T>();
        int_bar->_g_g_transpose = _g_g_transpose.template cast<T>();
        int_bar->_P = _P.template cast<T>();
        int_bar->use_low_rank_updates = use_low_rank_updates;

        return int_bar;
    };

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

    std::vector<std::vector<InterpolantDouble> > &get_basis() {
        return _unisolvent_basis;
    }

    Matrix get_P() {
        return _P;
    }

    void configure(pt::ptree & config);

    unsigned _max_polynomial_degree;
    unsigned _num_variable_symbols;
    //Turn _unisolvent_basis into Matrix
    std::vector<std::vector<InterpolantDouble> > _unisolvent_basis;
    Matrix _intermediate_matrix;
    Matrix _preintermediate_matrix;
    CustomLLT<Matrix, Eigen::Lower> _intermediate_LLT;
    Matrix _Q;
    Matrix _V;
    unsigned _L, _U;

    //Weighted polynomials
    Vector _g;
    Matrix _g_g_transpose;
    Matrix _P;

    //Usage of Woodburry matrix identity and simple heuristic for choice of variables
    //to update.
    bool use_low_rank_updates = true;

private:


    void compute_V_transpose_V();

    void construct_univariate(Vector poly_g);
    void construct_bivariate(Vector poly_g);
    void construct_multivariate(Vector poly_g);


};

#include "InterpolantDualSOSBarrier.hpp"

#endif //SOS_INTERPOLANTDUALSOSBARRIER_H
