// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "InterpolantDualSOSBarrier.h"
#include <boost/math/special_functions/binomial.hpp>

InterpolantDualSOSBarrier::InterpolantDualSOSBarrier(unsigned max_polynomial_degree_, Vector poly_g, unsigned num_variable_symbols_)
        : _max_polynomial_degree(max_polynomial_degree_), _num_variable_symbols(num_variable_symbols_) {

    //TODO: Check if still true for multivariate case.
    assert(poly_g.rows() <= max_polynomial_degree_ + 1);

    //poly_g.rows() is degree + 1 of the polynomial g;
    //TODO: setting _L using the number of rows only works for univariate polynomials
    //
    
    //TODO: check if type cast is safe.
    _L =  static_cast<unsigned>(boost::math::binomial_coefficient<double>(
            _num_variable_symbols + _max_polynomial_degree  + 1 - (unsigned) poly_g.rows(), _num_variable_symbols));
    _U =  static_cast<unsigned>(boost::math::binomial_coefficient<double>(
            2 * _max_polynomial_degree + _num_variable_symbols, _num_variable_symbols));

    _preintermediate_matrix = Matrix(_U, _L);
    _intermediate_matrix = Matrix(_L, _L);
    _intermediate_LLT = Eigen::LLT<Matrix>(_L);
    _V = Matrix(_L, _U);
    _Q = Matrix(_V.cols(), _V.cols());

    _num_variables = _U;
    _unisolvent_basis.resize(_U);

    if(_num_variable_symbols == 1) {
        construct_univariate(poly_g);
    } else if (_num_variable_symbols == 2){
        construct_bivariate(poly_g);
    } else {
        construct_multivariate(poly_g);
    }

};

void InterpolantDualSOSBarrier::construct_univariate(Vector poly_g)
{
    for (unsigned i = 0; i < _unisolvent_basis.size(); ++i) {
        BoostDouble cos_i = boost::multiprecision::cos(i * boost::math::constants::pi<BoostDouble>() / (_U - 1));
        InterpolantDouble dummy_ipm;
        InterpolantDouble cos_val = static_cast<InterpolantDouble>(cos_i);
        _unisolvent_basis[i] = cos_val;
    }

    //TODO: Figure out how choice of P could influence condition / stability of maps.

    //Use monomial standard basis to orthogonalize

    //_P is used in the Interior Point Method. Therefore we need to convert the multi-precision
    // floating-point into the IPM floating point precision

    //TODO: This interpolant Matrix only needs to be found once and can then be reused;
    spdlog::info("Construct interpolant point Matrix P...");

    //Alternative approach of finding _P via Chebyshev basis
    Eigen::MatrixXd cheb_P = ChebTools::u_matrix_library.get(_U - 1).block(0, 0, _U, _L);

    _g = Vector::Zero(_U);
    for (int p = 0; p < _U; ++p) {
        _g(p) = poly_g(0);
        for (int i = 1; i < poly_g.rows(); i++) {
            _g(p) += poly_g(i) * pow(_unisolvent_basis[p], i).convert_to<IPMDouble>();
        }
    }
    _g_g_transpose = _g * _g.transpose();

    //TODO: Option to compute cheb_P via InterpolantDouble;

    //TODO: Figure out whether orthogonalization could be done in double precision to speed up initialisation.
    std::cout << "Constructed." << std::endl;
    std::cout << "Orthogonalize..." << std::endl;
    cxxtimer::Timer orth_timer;
    orth_timer.start();
    Matrix P_tmp = cheb_P.cast<IPMDouble>();
    Matrix P_ortho = P_tmp.householderQr().householderQ();
    P_ortho.colwise().hnormalized();
    _P = P_ortho.block(0, 0, _U, _L);

    orth_timer.stop();
    std::cout << "Orthogonalization done in " << orth_timer.count<std::chrono::milliseconds>() / 1000.
              << " seconds." << std::endl;
}

void InterpolantDualSOSBarrier::construct_bivariate(Vector poly_g){

}

void InterpolantDualSOSBarrier::construct_multivariate(Vector poly_g){

}


bool InterpolantDualSOSBarrier::update_gradient_hessian_LLT(Vector x, bool check_interior_only) {
    _preintermediate_matrix.noalias() = _g.cwiseProduct(x).asDiagonal() * _P;
    _intermediate_matrix.noalias() = _P.transpose() * _preintermediate_matrix;
    _intermediate_LLT.compute(_intermediate_matrix);

    if (_intermediate_LLT.info() == Eigen::NumericalIssue) {
        return false;
    }

    if (check_interior_only) {
        return true;
    }

    _V.noalias() = _intermediate_LLT.matrixL().solve(_P.transpose());
    //Experiments showed that using the triangularView instead would slow down the program.
    //So we use the full Matrix _V to compute its product with the transpose.
    _Q.noalias() = _V.transpose() * _V;

    //TODO: store hessian as self-adjoint
    if (_stored_hessians.empty()) {
        _stored_hessians.resize(1);
    }
    _stored_hessians[0].first = x;
    _stored_hessians[0].second.noalias() = _g_g_transpose.cwiseProduct(_Q.cwiseProduct(_Q));

    if (_stored_gradients.empty()) {
        _stored_gradients.resize(1);
    }
    _stored_gradients[0].first = x;
    _stored_gradients[0].second.noalias() = -_Q.diagonal().cwiseProduct(_g);

    if (_stored_LLT.empty()) {
        _stored_LLT.resize(1);
    }

    _stored_LLT[0].first = x;
    _stored_LLT[0].second.compute(_stored_hessians[0].second.selfadjointView<Eigen::Lower>());

    return true;
}

Vector InterpolantDualSOSBarrier::gradient(Vector x) {
    auto *grad_ptr = find_gradient(x);
    if (grad_ptr) {
        return *grad_ptr;
    }
    update_gradient_hessian_LLT(x);
    return _stored_gradients[0].second;
}

Matrix InterpolantDualSOSBarrier::hessian(Vector x) {
    auto *hess_ptr = find_hessian(x);
    if (hess_ptr) {
        return *hess_ptr;
    }
    update_gradient_hessian_LLT(x);
    return _stored_hessians[0].second;
}

Eigen::LLT<Matrix> InterpolantDualSOSBarrier::llt(Vector x, bool) {
    auto *llt_ptr = find_LLT(x);
    if (llt_ptr) {
        return *llt_ptr;
    }
    update_gradient_hessian_LLT(x);
    return _stored_LLT[0].second;
}

//Should not be invoked as it is slow.
Matrix InterpolantDualSOSBarrier::inverse_hessian(Vector x) {
    Matrix L_inv = llt(x).matrixL().toDenseMatrix().inverse();
    Matrix inv = L_inv.transpose() * L_inv;
    return inv;
}

bool InterpolantDualSOSBarrier::in_interior(Vector x) {
    //The computational effort to calculate whether x is in the interior is nearly as high
    //as computing gradient and hessian. Therefore we just calculate them here as well.
    bool check_interior_only = false;
    return update_gradient_hessian_LLT(x, check_interior_only);
}

IPMDouble InterpolantDualSOSBarrier::concordance_parameter(Vector) {
    return _L;
}

Vector InterpolantDualSOSBarrier::initialize_x() {
    return Vector::Ones(_U);
}

Vector InterpolantDualSOSBarrier::initialize_s() {
    return -gradient(initialize_x());
}

