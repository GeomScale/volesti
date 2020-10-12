// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "InterpolantDualSOSBarrier.h"
#include <boost/math/special_functions/binomial.hpp>
#include "Padua/padua.h"

template <typename IPMDouble>
InterpolantDualSOSBarrier<IPMDouble>::InterpolantDualSOSBarrier(
        unsigned max_polynomial_degree_, Vector poly_g,unsigned num_variable_symbols_)
        : _max_polynomial_degree(max_polynomial_degree_), _num_variable_symbols(num_variable_symbols_) {

    //TODO: Check if still true for multivariate case.
    assert(poly_g.rows() <= max_polynomial_degree_ + 1);

    this->_custom_timers.resize(10);

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
    _intermediate_LLT = CustomLLT<Matrix, Eigen::Lower>(_L);
    _V = Matrix(_L, _U);
    _Q = Matrix(_V.cols(), _V.cols());

    this->_num_variables = _U;
    _unisolvent_basis.resize(_U);

    if (_num_variable_symbols == 1) {
        construct_univariate(poly_g);
    } else if (_num_variable_symbols == 2) {
        construct_bivariate(poly_g);
    } else {
        construct_multivariate(poly_g);
    }

};

template<typename IPMDouble>
void InterpolantDualSOSBarrier<IPMDouble>::construct_univariate(Vector poly_g) {
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

    //Alternative approach of finding _P via Chebyshev basis

    //TODO:Make this library more precise

    //Computing _g should be done with transformation matrix.
    _g = Vector::Zero(_U);
    for (int p = 0; p < _U; ++p) {
        _g(p) = poly_g(0);
        for (int i = 1; i < poly_g.rows(); i++) {
            _g(p) += poly_g(i) * pow(_unisolvent_basis[p][0], i).template convert_to<IPMDouble>();
        }
    }
    _g_g_transpose = _g * _g.transpose();

    //TODO: Option to compute cheb_P via InterpolantDouble;

    //TODO: Figure out whether orthogonalization could be done in double precision to speed up initialisation.
    spdlog::info("Construct orthogonal interpolant point Matrix P...");
    cxxtimer::Timer orth_timer;
    orth_timer.start();
    _P = orthogonal_P_Matrix_library<IPMDouble>.get(_L,_U);
    orth_timer.stop();
    std::cout << "Orthogonalization done in " << orth_timer.count<std::chrono::milliseconds>() / 1000.
              << " seconds." << std::endl;
}

//Untested
template<typename IPMDouble>
void InterpolantDualSOSBarrier<IPMDouble>::construct_bivariate(Vector poly_g) {

    //For now no weighted polynomials
    assert(poly_g == Vector::Ones(1));
    unsigned const corrected_d = 2 * _max_polynomial_degree + 1 - poly_g.size();
    unsigned const corrected_d_plus_1 = corrected_d + 1;

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
template<typename IPMDouble>
void InterpolantDualSOSBarrier<IPMDouble>::construct_multivariate(Vector poly_g) {
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

    DegreeTuple dt(_num_variable_symbols, 2 * _max_polynomial_degree);
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
template<typename IPMDouble>
void InterpolantDualSOSBarrier<IPMDouble>::compute_V_transpose_V() {
    _Q.noalias() = _V.transpose() * _V;
}

template<typename IPMDouble>
void InterpolantDualSOSBarrier<IPMDouble>::configure(pt::ptree & config){
    if(config.find("use_low_rank_updates") != config.not_found()){
        use_low_rank_updates = config.get<bool>("use_low_rank_updates");
    }
}

template<typename IPMDouble>
bool InterpolantDualSOSBarrier<IPMDouble>::update_gradient_hessian_LLT(Vector x, bool check_interior_only) {

    Matrix Q;
    CustomLLT<Matrix, Eigen::Lower> LLT;
    IPMDouble stored_gx_norm;
    Vector new_stored_x;

//    if(!this->_stored_gradients.empty()){
//        //This is to test whether lower ran updates make sense
//
//        Vector stored_x = this->_stored_gradients[0].first;
//        stored_gx_norm = _g.cwiseProduct(stored_x).norm();
//        Matrix aux(2, x.rows());
//        IPMDouble relative_norm = stored_x.norm()/x.norm();
////        std::cout << "relative norm is " << relative_norm << std::endl;
//        aux << stored_x.transpose(), x.transpose() * relative_norm;
////        std::cout << "aux is\n" << aux << std::endl;
//        Vector rel = stored_x.cwiseProduct(x.cwiseInverse());
//        std::vector<IPMDouble> relative_vec(rel.data(), rel.data() + rel.rows());
//        sort(relative_vec.begin(), relative_vec.end());
//        IPMDouble mean = relative_vec[relative_vec.size() / 2];
//
//        std::cout << "Relative mean vector is: \n";
//
//        unsigned small_variation_count_05 = 0;
//        unsigned small_variation_count_01 = 0;
//        for(IPMDouble i : relative_vec){
//            i/=mean;
//            std::cout << i << ", ";
//            if(abs(i-1) < .05){
//                small_variation_count_05++;
//            }
//            if(abs(i-1) < .01){
//                small_variation_count_01++;
//            }
//        }
//        std::cout << std::endl;
//        std::cout << "small variation ratio: " << double(small_variation_count_01) / relative_vec.size() << std::endl;
//        std::cout << "small variation ratio: " << double(small_variation_count_05) / relative_vec.size() << std::endl;
//    }

    bool do_exact_computation = true;

    if (use_low_rank_updates and not this->_stored_gradients.empty()) {
        do_exact_computation = false;
        this->_custom_timers[9].start();
        Vector stored_x = this->_stored_gradients[0].first;
        stored_gx_norm = _g.cwiseProduct(stored_x).norm();
        Matrix aux(2, x.rows());
        aux << stored_x.transpose(), x.transpose();
        std::cout << "aux is\n" << aux << std::endl;
        Vector stored_scaled_gx = _g.cwiseProduct(stored_x).normalized();
        //TODO: Find best scaling to minimize number of adjusted variables.
        Vector scaled_gx = _g.cwiseProduct(x).normalized();
        Vector relative = scaled_gx - stored_scaled_gx;
        Vector relative_abs = relative.cwiseAbs();
//        TODO: Check if adding diagonal here makes sense. This corresponds to the gradient and the term occurs for the Sherman-Morrison update.

        std::vector<IPMDouble> relative_vec(relative_abs.data(), relative_abs.data() + relative_abs.rows());
        std::vector<size_t> relative_vec_sorted = sort_indexes(relative_vec);
        std::reverse(relative_vec_sorted.begin(), relative_vec_sorted.end());

        Q = _Q;
        LLT = _intermediate_LLT;
        assert(_U == relative_vec_sorted.size());

        new_stored_x = this->_stored_gradients[0].first;
        for (int i = 0; i < _U; i++) {

            //rank one updates

            unsigned idx = relative_vec_sorted[i];

            //TODO: find good threshold
            if (relative_abs(idx) < 1e-10 * relative_abs(relative_vec_sorted[0])) {
                std::cout << "break after " << i << " of " << _U << " updates"
                << std::endl;
                break;
            }

            //Q update
            new_stored_x(idx) = x(idx);

            this->_custom_timers[7].start();
            Vector P_idx = _P.row(idx).transpose();
            Vector tmp = LLT.solve(P_idx);
            Vector rank_one_factor = _P * tmp;
            Matrix Q_update =
                    -rank_one_factor * rank_one_factor.transpose() /
                    (1 / (stored_gx_norm * relative(idx)) + tmp.dot(P_idx));
//            std::cout << "Q_update normed is \n" << Q_update  << std::endl;

            Q += Q_update;
            this->_custom_timers[7].stop();

            //L update
            this->_custom_timers[8].start();
            LLT.rankUpdate(P_idx, stored_gx_norm * relative(idx));
            if (LLT.info() == Eigen::NumericalIssue) {
                //This means we should either change the order in which we apply the rank-1 updates
                //or revert back to exact computation. Currently we just jump to exact computation.

                do_exact_computation = true;
                break;
            }
            this->_custom_timers[8].stop();
        }
        this->_custom_timers[9].stop();

        if (not do_exact_computation) {
//        std::cout << "Old matrix LLT \n" << LLT.matrixL().toDenseMatrix() << std::endl;
            _intermediate_LLT.copy_and_scale(LLT, sqrt(_g.cwiseProduct(x).norm() / stored_gx_norm));
//        std::cout << "Newly copied matrix \n"
//        << _intermediate_LLT.matrixL().toDenseMatrix() / sqrt(_g.cwiseProduct(x).norm() / stored_gx_norm) << std::endl;
            _Q = Q * stored_gx_norm / _g.cwiseProduct(x).norm();
        }
    }

    if (do_exact_computation) {
        new_stored_x = x;

        if(this->_stored_intermediate_LLT.empty()){
            this->_stored_intermediate_LLT.resize(1);
            this->_stored_intermediate_LLT[0].first = Vector::Zero(new_stored_x.rows());
        }

        if(this->_stored_intermediate_LLT[0].first == new_stored_x){
            _intermediate_LLT = this->_stored_intermediate_LLT[0].second;
        } else {

            this->_custom_timers[1].start();
            _preintermediate_matrix.noalias() = _g.cwiseProduct(x).asDiagonal() * _P;
            _intermediate_matrix.noalias() = _P.transpose() * _preintermediate_matrix;
            this->_custom_timers[1].stop();
            this->_custom_timers[2].start();
            _intermediate_LLT.compute(_intermediate_matrix);
            this->_custom_timers[2].stop();

            if (_intermediate_LLT.info() == Eigen::NumericalIssue) {
                return false;
            } else {
                this->_stored_intermediate_LLT[0].first = new_stored_x;
                this->_stored_intermediate_LLT[0].second = _intermediate_LLT;
            }
        }

        if (check_interior_only) {
            return true;
        }

        this->_custom_timers[3].start();
        _V.noalias() = _intermediate_LLT.matrixL().solve(_P.transpose());
        this->_custom_timers[3].stop();

        //Experiments showed that using the triangularView instead would slow down the program.
        //So we use the full Matrix _V to compute its product with the transpose.

        this->_custom_timers[4].start();
        compute_V_transpose_V();
        this->_custom_timers[4].stop();
    }

    //TODO: store hessian as self-adjoint
    if (this->_stored_hessians.empty()) {
        this->_stored_hessians.resize(1);
    }
    this->_stored_hessians[0].first = new_stored_x;
    this->_stored_hessians[0].second.noalias() = _g_g_transpose.cwiseProduct(_Q.cwiseProduct(_Q));

    if (this->_stored_gradients.empty()) {
        this->_stored_gradients.resize(1);
    }
    this->_stored_gradients[0].first = new_stored_x;
    this->_stored_gradients[0].second.noalias() = -_Q.diagonal().cwiseProduct(_g);

    if (this->_stored_LLT.empty()) {
        this->_stored_LLT.resize(1);
    }

    this->_stored_LLT[0].first = new_stored_x;
    this->_custom_timers[5].start();
    this->_stored_LLT[0].second.compute(this->_stored_hessians[0].second.template selfadjointView<Eigen::Lower>());
    this->_custom_timers[5].stop();

    return true;
}

template<typename IPMDouble>
Vector<IPMDouble> InterpolantDualSOSBarrier<IPMDouble>::gradient(Vector x) {
    auto *grad_ptr = this->find_gradient(x);
    if (grad_ptr) {
        return *grad_ptr;
    }
    update_gradient_hessian_LLT(x);
    return this->_stored_gradients[0].second;
}

template<typename IPMDouble>
Matrix<IPMDouble> InterpolantDualSOSBarrier<IPMDouble>::hessian(Vector x) {
    auto *hess_ptr = this->find_hessian(x);
    if (hess_ptr) {
        return *hess_ptr;
    }
    update_gradient_hessian_LLT(x);
    return this->_stored_hessians[0].second;
}

template<typename IPMDouble>
Eigen::LLT<Matrix<IPMDouble> > InterpolantDualSOSBarrier<IPMDouble>::llt(Vector x, bool) {
    auto *llt_ptr = this->find_LLT(x);
    if (llt_ptr) {
        return *llt_ptr;
    }
    update_gradient_hessian_LLT(x);
    return this->_stored_LLT[0].second;
}

//Should not be invoked as it is slow.
template<typename IPMDouble>
Matrix<IPMDouble> InterpolantDualSOSBarrier<IPMDouble>::inverse_hessian(Vector x) {
    Matrix L_inv = llt(x).matrixL().toDenseMatrix().inverse();
    Matrix inv = L_inv.transpose() * L_inv;
    return inv;
}

template<typename IPMDouble>
bool InterpolantDualSOSBarrier<IPMDouble>::in_interior(Vector x) {
    //The computational effort to calculate whether x is in the interior is nearly as high
    //as computing gradient and hessian. Therefore we just calculate them here as well.
    bool check_interior_only = true;
    return update_gradient_hessian_LLT(x, check_interior_only);
}

template<typename IPMDouble>
IPMDouble InterpolantDualSOSBarrier<IPMDouble>::concordance_parameter(Vector) {
    return _L;
}

template<typename IPMDouble>
Vector<IPMDouble> InterpolantDualSOSBarrier<IPMDouble>::initialize_x() {
    return Vector::Ones(_U);
}

template<typename IPMDouble>
Vector<IPMDouble> InterpolantDualSOSBarrier<IPMDouble>::initialize_s() {
    return -gradient(initialize_x());
}

