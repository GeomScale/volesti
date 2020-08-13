// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "InterpolantDualSOSBarrier.h"
#include <boost/math/special_functions/binomial.hpp>
#include "../../../external/Padua/padua.h"

InterpolantDualSOSBarrier::InterpolantDualSOSBarrier(unsigned max_polynomial_degree_, Vector poly_g,
                                                     unsigned num_variable_symbols_)
        : _max_polynomial_degree(max_polynomial_degree_), _num_variable_symbols(num_variable_symbols_) {

    //TODO: Check if still true for multivariate case.
    assert(poly_g.rows() <= max_polynomial_degree_ + 1);

    //poly_g.rows() is degree + 1 of the polynomial g;
    //TODO: setting _L using the number of rows only works for univariate polynomials
    //

    //For now no weights for multivariate case.
    assert(_num_variable_symbols == 1 or poly_g == Vector::Ones(poly_g.rows()));

    //TODO: check if type cast is safe.
    _L = static_cast<unsigned>(boost::math::binomial_coefficient<double>(
            _num_variable_symbols + _max_polynomial_degree + 1 - (unsigned) poly_g.rows(), _num_variable_symbols));
    _U = static_cast<unsigned>(boost::math::binomial_coefficient<double>(
            2 * _max_polynomial_degree + _num_variable_symbols, _num_variable_symbols));

    _preintermediate_matrix = Matrix(_U, _L);
    _intermediate_matrix = Matrix(_L, _L);
    _intermediate_LLT = Eigen::LLT<Matrix>(_L);
    _V = Matrix(_L, _U);
    _Q = Matrix(_V.cols(), _V.cols());

    _num_variables = _U;
    _unisolvent_basis.resize(_U);

    if (_num_variable_symbols == 1) {
        construct_univariate(poly_g);
    } else if (_num_variable_symbols == 2) {
        construct_bivariate(poly_g);
    } else {
        construct_multivariate(poly_g);
    }

};

void InterpolantDualSOSBarrier::construct_univariate(Vector poly_g) {
    for (unsigned i = 0; i < _unisolvent_basis.size(); ++i) {
        BoostDouble cos_i = boost::multiprecision::cos(i * boost::math::constants::pi<BoostDouble>() / (_U - 1));
        InterpolantDouble dummy_ipm;
        InterpolantDouble cos_val = static_cast<InterpolantDouble>(cos_i);
        _unisolvent_basis[i].push_back(cos_val);
    }

    //TODO: Figure out how choice of P could influence condition / stability of maps.

    //Use monomial standard basis to orthogonalize

    //_P is used in the Interior Point Method. Therefore we need to convert the multi-precision
    // floating-point into the IPM floating point precision

    //TODO: This interpolant Matrix only needs to be found once and can then be reused;
    spdlog::info("Construct interpolant point Matrix P...");

    //Alternative approach of finding _P via Chebyshev basis

    //TODO:Make this library more precise
    Eigen::MatrixXd cheb_P = ChebTools::u_matrix_library.get(_U - 1).block(0, 0, _U, _L);

    //Computing _g should be done with transformation matrix.
    _g = Vector::Zero(_U);
    for (int p = 0; p < _U; ++p) {
        _g(p) = poly_g(0);
        for (int i = 1; i < poly_g.rows(); i++) {
            _g(p) += poly_g(i) * pow(_unisolvent_basis[p][0], i).convert_to<IPMDouble>();
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

//Untested
void InterpolantDualSOSBarrier::construct_bivariate(Vector poly_g) {

    //Error waiting to happen. Here the unisolvent basis is correctly used for _L, but in the univariate
    //case we use it for a basis or size _U.
//    _unisolvent_basis.resize(_L);
    //For now no weighted polynomials
    assert(poly_g == Vector::Ones(1));
    unsigned const corrected_d = 2 * _max_polynomial_degree + 1 - poly_g.size();
    unsigned const corrected_d_plus_1 = corrected_d + 1;

    //Order of bivariate unisolvent basis: first points according to Even x Odd, then Odd x Even.
//    std::vector<InterpolantDouble> cos_values_d;
//    std::vector<InterpolantDouble> cos_values_d_plus_1;
//
//    for (unsigned i = 0; i <= corrected_d; i++) {
//        BoostDouble cos_i = boost::multiprecision::cos(i * boost::math::constants::pi<BoostDouble>() / corrected_d);
//        InterpolantDouble cos_val = static_cast<InterpolantDouble>(cos_i);
//        cos_values_d.push_back(cos_val);
//    }
//
//    for (unsigned i = 0; i <= corrected_d_plus_1; i++) {
//        BoostDouble cos_i = boost::multiprecision::cos(i * boost::math::constants::pi<BoostDouble>() / corrected_d_plus_1);
//        InterpolantDouble cos_val = static_cast<InterpolantDouble>(cos_i);
//        cos_values_d_plus_1.push_back(cos_val);
//    }
//
//    unsigned idx = 0;
//    for (int i = 0; i <= corrected_d; i+=2) {
//       for(int j = 1; j <= corrected_d_plus_1; j+=2){
//           _unisolvent_basis[idx++] = {cos_values_d[i], cos_values_d_plus_1[j]};
//       }
//    }
//
//    for (int i = 1; i <= corrected_d; i+=2) {
//        for(int j = 0; j <= corrected_d_plus_1; j+=2){
//            _unisolvent_basis[idx++] = {cos_values_d[i], cos_values_d_plus_1[j]};
//        }
//    }

    double *pd_pts = padua::padua_points(_max_polynomial_degree + 1);

    for (int i = 0; i < _U; i++) {
        _unisolvent_basis[i] = {pd_pts[2 * i], pd_pts[2 * i + 1]};
    }

    //Set weight vector _g;
    //TODO: do properly for weighted case.
    _g = Vector::Ones(_U);
    _g_g_transpose = _g * _g.transpose();

    std::cout << "Construct Matrix P" << std::endl;

    _P.resize(_U, _L);

    //TODO: Make following loops more efficient OR find mathematical theory that simplifies expressions.
    unsigned col_idx = 0;
    for (int i = 0; i <= _max_polynomial_degree; i++) {
        for (int j = 0; j + i <= _max_polynomial_degree; j++) {
            Eigen::VectorXd vec_i = Eigen::VectorXd::Zero(i + 1);
            vec_i(i) = 1;
            Eigen::VectorXd vec_j = Eigen::VectorXd::Zero(j + 1);
            vec_j(j) = 1;
            ChebTools::ChebyshevExpansion cheb_i(vec_i);
            ChebTools::ChebyshevExpansion cheb_j(vec_j);
            for (int k = 0; k < _U; k++) {
                double first_eval = cheb_i.y_recurrence(static_cast<double>(_unisolvent_basis[k][0]));
                double second_eval = cheb_j.y_recurrence(static_cast<double>(_unisolvent_basis[k][1]));
                _P(k, col_idx) = first_eval * second_eval;
            }
            col_idx++;
        }
    }

    std::cout << "P before orthogonolisation: " << std::endl << _P << std::endl;

    Matrix P_ortho = _P.householderQr().householderQ();
    P_ortho.colwise().hnormalized();
    _P = P_ortho.block(0, 0, _U, _L);
    //end
}

//TODO: test
void InterpolantDualSOSBarrier::construct_multivariate(Vector poly_g) {
    assert(poly_g == Vector::Ones(1));

    //Set weight vector _g;
    //TODO: do properly for weighted case.
    _g = Vector::Ones(_U);
    _g_g_transpose = _g * _g.transpose();

    _P.resize(_U, _L);

    //just for testing.

    int num_candidates = 1;
    for (int i = 2 * _max_polynomial_degree + 1 + 1;
         i <= 2 * _max_polynomial_degree + _num_variable_symbols + 1; i++) {
        num_candidates *= i;
    }

    //generate Fekete candidates


    std::vector<std::vector<double> > candidates;
    std::vector<unsigned> comb_bound;
    for (int i = 2 * _max_polynomial_degree + 1; i <= 2 * _max_polynomial_degree + _num_variable_symbols; i++) {
        comb_bound.push_back(i);
    }
    AllCombinationTuple comb_tuple(comb_bound);

    do {
        std::vector<unsigned> &cheb_v = comb_tuple.get_combination();
        std::vector<double> cand;
        for (int i = 0; i < cheb_v.size(); i++) {
            cand.push_back(cos((double) cheb_v[i] * boost::math::constants::pi<double>() / (double) comb_bound[i]));
        }
        candidates.push_back(cand);
    } while (comb_tuple.next());

    for (auto cand : candidates) {
        for (auto c : cand) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }

    assert(num_candidates == candidates.size());

    //generate Fekete polynomials

    Matrix candidate_matrix(_U, candidates.size());

    DegreeTuple dt(_num_variable_symbols, 2*_max_polynomial_degree);
    unsigned tup_idx = 0;
    do {
        //Compute chebyshev polynomial evaluation
        std::vector<unsigned> &tup = dt.get_tuple();
        for (int i = 0; i < candidates.size(); i++) {
            double cheb_eval = 1.;
            for (int j = 0; j < candidates[i].size(); j++) {
                Eigen::VectorXd vec_j = Eigen::VectorXd::Zero(tup[j] + 1);
                vec_j(tup[j]) = 1;
                ChebTools::ChebyshevExpansion cheb_j(vec_j);
                double eval_j = cheb_j.y_recurrence(static_cast<double>(candidates[i][j]));
                cheb_eval *= eval_j;
            }
            candidate_matrix(tup_idx, i) = cheb_eval;
        }
        tup_idx++;
    } while (dt.next_valid());


    assert(tup_idx == _U);

//    std::cout << "Candidate matrix is \n";
//    std::cout << candidate_matrix << std::endl;


    Matrix col_permutated_cand_matrix = candidate_matrix * candidate_matrix.colPivHouseholderQr().colsPermutation();
//    std::cout << "And after permutation \n";
//    std::cout << col_permutated_cand_matrix << std::endl;
    _P = col_permutated_cand_matrix.block(0, 0, _U, _U).transpose();
}

//For profiling purposes
void InterpolantDualSOSBarrier::compute_V_transpose_V(){
        _Q.noalias() = _V.transpose() * _V;
        double test = _Q(1,1);
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

    compute_V_transpose_V();

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

