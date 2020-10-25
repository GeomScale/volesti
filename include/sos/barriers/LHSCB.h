// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SOS_LHSCB_H
#define SOS_LHSCB_H

#include<iostream>
#include <vector>
#include "utils.h"

template <typename IPMDouble>
class LHSCB {
    typedef Vector<IPMDouble> Vector;
    typedef Matrix<IPMDouble> Matrix;
public:
    LHSCB() : _num_variables(0) {};

    template <typename T>
    void cast_members_from(const LHSCB<T> & other){
       _num_variables = other._num_variables;
    }

    template<typename T>
    LHSCB<T> * cast(){
        LHSCB<T> * lhscb;
        (*lhscb).template cast_members_from<IPMDouble>(*this);
        return lhscb;
    }

    virtual ~LHSCB() {};

    virtual Vector gradient(Vector x) = 0;

    virtual Matrix hessian(Vector x) = 0;

    virtual Eigen::LLT<Matrix> llt(Vector x, bool symmetrize = 0);

    virtual Matrix llt_solve(Vector x, const Matrix &rhs);

    virtual Vector llt_L_solve(Vector x, Vector rhs);

    Vector *find_gradient(Vector x);

    Matrix *find_hessian(Vector x);

    Eigen::LLT<Matrix> *find_LLT(Vector x);

    virtual Matrix inverse_hessian(Vector x);

    virtual bool in_interior(Vector x) = 0;

    virtual IPMDouble concordance_parameter(Vector x) = 0;

    virtual Vector initialize_x(IPMDouble parameter) {
        return parameter * initialize_x();
    }

    //TODO: figure out if initializing the dual is in general just - 1/mu * g(x) (i.e. whether this is in the dual cone)
    virtual Vector initialize_s(IPMDouble parameter) {
        return initialize_s() / parameter;
    }

    virtual Vector initialize_x() = 0;

    virtual Vector initialize_s() = 0;

    cxxtimer::Timer _in_interior_timer;
    std::vector<cxxtimer::Timer> _custom_timers;

    unsigned _num_variables;


protected:

    std::vector<std::pair<Vector, Vector> > _stored_gradients;
    std::vector<std::pair<Vector, Matrix> > _stored_hessians;
    std::vector<std::pair<Vector, Eigen::LLT<Matrix> > > _stored_LLT;
    std::vector<std::pair<Vector, CustomLLT<Matrix, Eigen::Lower> > > _stored_intermediate_LLT;

public:
    unsigned int getNumVariables() const;
};

#include "LHSCB.hpp"

#endif //SOS_LHSCB_H
