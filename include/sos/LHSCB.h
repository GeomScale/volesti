// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_LHSCB_H
#define NONSYMMETRICCONICOPTIMIZATION_LHSCB_H

#include<iostream>
#include <vector>
#include "utils.cpp"
#include "ChebTools/ChebTools.h"



class LHSCB {
public:
    LHSCB() : _num_variables(0) {};

    virtual ~LHSCB() {};

    virtual Vector gradient(Vector x) = 0;

    virtual Matrix hessian(Vector x) = 0;

    virtual Eigen::LLT<Matrix> llt(Vector x, bool symmetrize = 0);

    virtual Matrix llt_solve(Vector x, const Matrix & rhs);

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


protected:
    unsigned _num_variables;
    std::vector<std::pair<Vector, Vector> > _stored_gradients;
    std::vector<std::pair<Vector, Matrix> > _stored_hessians;
    std::vector<std::pair<Vector, Eigen::LLT<Matrix> > > _stored_LLT;

public:
    unsigned int getNumVariables() const;
};

//Not sure about inverse Hessian. Might treat free variables differently.
//Should FullSpaceBarrier and ZeroSpaceBarrier even exist or should the corresponding variables be handled differently in the algorithm itself?

class FullSpaceBarrier final : public LHSCB {

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    Matrix inverse_hessian(Vector x) override;

    bool in_interior(Vector x) override;

    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

public:
    FullSpaceBarrier(unsigned num_variables_) {
        _num_variables = num_variables_;
    }
};

class ZeroSpaceBarrier final : public LHSCB {

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

    Matrix inverse_hessian(Vector x) override;

    bool in_interior(Vector x) override;

    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

public:
    ZeroSpaceBarrier(unsigned num_variables_) {
        _num_variables = num_variables_;
    }
};

class LPStandardBarrier final : public LHSCB {
public:
    LPStandardBarrier() : LHSCB() {};

    LPStandardBarrier(unsigned num_variables_) {
        _num_variables = num_variables_;
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

//FIXME: Class currently only supports univariate problems
//TODO: set variable _num_variables at different spot. Currently it is initialized in each individual constructor.
// Write method that gets invoked by parent constructor.

//TODO: Add WSOS barrier.

class InterpolantDualSOSBarrier : public LHSCB {

public:
    InterpolantDualSOSBarrier() : LHSCB() {};

    InterpolantDualSOSBarrier(unsigned max_polynomial_degree_) : _max_polynomial_degree(max_polynomial_degree_) {

        _L = _max_polynomial_degree + 1;
        _U = 2 * _max_polynomial_degree + 1;
        _num_variables = _U;
        _unisolvent_basis.resize(_U);

        for (unsigned i = 0; i < _unisolvent_basis.size(); ++i) {
            BoostDouble cos_i = boost::multiprecision::cos(i * boost::math::constants::pi<BoostDouble>() / (_U - 1));
            InterpolantDouble dummy_ipm;
            InterpolantDouble cos_val = InterpolantDoubletoIPMDouble(cos_i, dummy_ipm);
            _unisolvent_basis[i] = cos_val;
        }

//        InterpolantMatrix P_interp = InterpolantMatrix(_U, _L);

        //TODO: Figure out how choice of P could influence condition / stability of maps.

        //Use monomial standard basis to orthogonalize

        //_P is used in the Interior Point Method. Therefore we need to convert the multi-precision
        // floating-point into the IPM floating point precision

        //TODO: This interpolant Matrix only needs to be found once and can then be reused;
        std::cout << "Construct interpolant point Matrix P..." << std::endl;

        //Construct P_interp via Chebyshev basis.

//        for (unsigned row = 0; row < P_interp.rows(); ++row) {
//            for (unsigned col = 0; col < P_interp.cols(); ++col) {
//                //TODO: make more efficient.
//                P_interp(row, col) = col == 0 ? InterpolantDouble(1.) : pow(_unisolvent_basis[row], col);
//            }
//        }

        //Alternative approach of finding _P via Chebyshev basis
        Eigen::MatrixXd cheb_P = ChebTools::u_matrix_library.get(_U - 1).block(0,0,_U, _L);

//        assert(cheb_P.cols() == _L and P_interp.cols() == _L);
//        Eigen::MatrixXd P_double = InterpolantMatrixToMatrix(P_interp, cheb_P);

//        std::cout << "P_double matrix is : \n" << P_double << std::endl;
//        std::cout << "cheb_P matrix is: \n" << cheb_P << std::endl;

//        Eigen::MatrixXd concat(P_double.rows(), 2 * _L);
//        concat << P_double, cheb_P;

//        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(concat);
//        auto rank = lu_decomp.rank();
//        std::cout << "_L is " << _L << ", _U is " <<  _U << ", and rank is " << rank << std::endl;


        //TODO: Option to compute cheb_P via InterpolantDouble;

        //TODO: Figure out whehter orthogonalization could be done in double precision to speed up initialisation.
        std::cout << "Constructed." << std::endl;
        std::cout << "Orthogonalize..." << std::endl;
        cxxtimer::Timer orth_timer;
        orth_timer.start();
//        if(not orthogonalize_in_ipm_double) {
//            InterpolantMatrix P_ortho = P_interp.householderQr().householderQ();
//            P_ortho.colwise().hnormalized();
//            InterpolantMatrix P_intermediate = P_ortho.block(0, 0, _U, _L);
//            _P = InterpolantMatrixToMatrix(P_intermediate, _P);
//        }
//        else {
//            Matrix P_tmp = InterpolantMatrixToMatrix(P_interp, _P);
            Matrix P_tmp = cheb_P.cast<IPMDouble>();
            Matrix P_ortho = P_tmp.householderQr().householderQ();
            P_ortho.colwise().hnormalized();
            _P = P_ortho.block(0,0,_U,_L);
//        }

//        std::cout << "Col product of cheb matrix: \n" << cheb_P.transpose() * cheb_P << std::endl;

        orth_timer.stop();
        std::cout << "Orthogonalization done in " << orth_timer.count<std::chrono::milliseconds>() / 1000.
                  << " seconds." << std::endl;

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

    std::vector<InterpolantDouble> &get_basis() {
        return _unisolvent_basis;
    }

private:
    unsigned _max_polynomial_degree;
    std::vector<InterpolantDouble> _unisolvent_basis;
    unsigned _L, _U;
    Matrix _P;

    //This variable sets whether the orthogonalisation of P should be done in either the Interpolant Double Type
    // or the IPM Basis Type
    bool orthogonalize_in_ipm_double = true;
};

class ProductBarrier : public LHSCB {

public:
    ProductBarrier() : LHSCB() {}

    ProductBarrier(std::vector<LHSCB *> barriers_, std::vector<unsigned> num_variables_) {
        assert(barriers_.size() == num_variables_.size());
        for (unsigned j = 0; j < barriers_.size(); ++j) {
            _barriers.push_back(barriers_[j]);
            _num_vars_per_barrier.push_back(num_variables_[j]);
        }
    }



    void add_barrier(LHSCB *lhscb) {
        _barriers.push_back(lhscb);
        _num_vars_per_barrier.push_back(lhscb->getNumVariables());
        _num_variables += lhscb->getNumVariables();

    }

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

//    Matrix inverse_hessian(Vector x) override;
    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

    Eigen::LLT<Matrix> llt(Vector x, bool symmetrize = 0) override;

    Matrix llt_solve(Vector x, const Matrix & rhs) override;

    Vector llt_L_solve(Vector x, Vector rhs) override;

    Matrix inverse_hessian(Vector x) override;

private:
    std::vector<LHSCB *> _barriers;
    std::vector<unsigned> _num_vars_per_barrier;
};

class DualSOSConeBarrier : public LHSCB {
    //Implementation for Standard Monomial Basis

public:
    DualSOSConeBarrier() : LHSCB() {};

    DualSOSConeBarrier(unsigned max_polynomial_degree_) : _max_polynomial_degree(max_polynomial_degree_) {
    };

    Vector gradient(Vector x) override;

    Matrix hessian(Vector x) override;

//    Matrix inverse_hessian(Vector x) override;
    bool in_interior(Vector x) override;

    //TODO: better solution for concordance parameter;
    IPMDouble concordance_parameter(Vector x) override;

    Vector initialize_x() override;

    Vector initialize_s() override;

    Matrix Lambda(Vector x);

private:
    unsigned _max_polynomial_degree;
};

#endif //NONSYMMETRICCONICOPTIMIZATION_LHSCB_H
