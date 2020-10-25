// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file
#ifndef SOS_NONSYMMETRICIPM_H
#define SOS_NONSYMMETRICIPM_H

#include "barriers/LHSCB.h"
#include "barriers/ProductBarrier.h"
#include <iostream>
#include <cxxtimer.hpp>
#include <fstream>


template<class T>
class Instance {
public:
    Constraints<T> constraints;
    LHSCB<T> *barrier;
};

template<class T>
class DirectionDecomposition {
public:
    DirectionDecomposition(Vector<T> v, unsigned const n, unsigned const m) {
        assert(v.rows() == m + n + 1 + n + 1);
        y = v.block(0, 0, m, 1);
        x = v.block(m, 0, n, 1);
        tau = v.block(m + n, 0, 1, 1).sum();
        s = v.block(m + n + 1, 0, n, 1);
        kappa = v.block(m + n + 1 + n, 0, 1, 1).sum();
    }

    friend std::ostream &operator<<(std::ostream &os, const DirectionDecomposition<T> &dir);

    Vector<T> y, x, s;
    T kappa, tau;
};

template<class T>
class ErrorConstants {
public:

    ErrorConstants() {};

    template<class U>
    ErrorConstants(const ErrorConstants<U> other) {
        //Note: Boost Dependency
        primal = boost::numeric_cast<T>(other.primal);
        dual = boost::numeric_cast<T>(other.dual);
        complementary = boost::numeric_cast<T>(other.complementary);
    }

    void set(Matrix<T> A, Vector<T> b, Vector<T> c) {
        //Constant are used to conform with the implementation by Papp & Yildiz.
        primal = std::max(T(1), sqrt(pow(A.norm(), 2) + pow(b.norm(), 2)));
        Matrix<T> Id = Matrix<T>::Identity(A.cols(), A.cols());
        dual = std::max(T(1), sqrt(pow(A.transpose().norm(), 2) + pow(Id.norm(), 2) + pow(c.norm(), 2)));
        complementary = sqrt(pow(c.norm(), 2) + pow(b.norm(), 2) + 1);
    }

    T primal;
    T dual;
    T complementary;
};

template<typename IPMDouble>
class NonSymmetricIPM {

    typedef Matrix<IPMDouble> Matrix;
    typedef Vector<IPMDouble> Vector;

public:

    NonSymmetricIPM(Matrix &A_, Vector &b_, Vector &c_, LHSCB<IPMDouble> *barrier_);

    NonSymmetricIPM(Instance<IPMDouble> &instance) : NonSymmetricIPM(instance.constraints.A,
                                                                     instance.constraints.b, instance.constraints.c,
                                                                     instance.barrier) {
        _logger->set_level(spdlog::level::info);
    };

    NonSymmetricIPM(Instance<IPMDouble> &, std::string);


    template<typename T>
    void cast_members_from(const NonSymmetricIPM<T> &other) {
        A_sparse = other.A_sparse.template cast<IPMDouble>();
        _basis_ker_A = other._basis_ker_A.template cast<IPMDouble>();
        _config = other._config;

        _num_predictor_steps = other._num_predictor_steps;
        _num_corrector_steps = other._num_corrector_steps;
        _param_step_length_predictor = boost::numeric_cast<IPMDouble>(other._param_step_length_predictor);
        _step_length_predictor = boost::numeric_cast<IPMDouble>(other._step_length_predictor);
        _step_length_corrector = boost::numeric_cast<IPMDouble>(other._step_length_corrector);
        _epsilon = boost::numeric_cast<T>(other._epsilon);

        _large_neighborhood = boost::numeric_cast<IPMDouble>(other._large_neighborhood);
        _small_neighborhood = boost::numeric_cast<IPMDouble>(other._small_neighborhood);
        //copy current solution

        x = other.x.template cast<IPMDouble>();
        y = other.y.template cast<IPMDouble>();
        s = other.s.template cast<IPMDouble>();
        kappa = boost::numeric_cast<IPMDouble>(other.kappa);
        tau = boost::numeric_cast<IPMDouble>(other.tau);

        _last_predictor_direction = other._last_predictor_direction.template cast<IPMDouble>();

        _err_consts = ErrorConstants<IPMDouble>(other._err_consts);

        //copy configuration;

        _use_line_search = other._use_line_search;
    }

    template<typename T>
    NonSymmetricIPM<T> *cast_with_product_barrier() {

        //cast initialisatino data
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A_ = A.template cast<T>();
        Eigen::Matrix<T, Eigen::Dynamic, 1> b_ = b.template cast<T>();
        Eigen::Matrix<T, Eigen::Dynamic, 1> c_ = c.template cast<T>();

        //TODO: High priority. Figure out how to undo the cast to the ProductBarrier.
        ProductBarrier<T> *barrier_ = static_cast<ProductBarrier<IPMDouble> *>(_barrier)->template cast<T>();

        assert(barrier_ != nullptr);

        NonSymmetricIPM<T> *nonSymmetricIPM = new NonSymmetricIPM<T>(A_, b_, c_, barrier_);

        //cast members

        nonSymmetricIPM->template cast_members_from<IPMDouble>(*this);

        //copy configuration;

        assert(nonSymmetricIPM->x.rows() == nonSymmetricIPM->s.rows());
        assert(nonSymmetricIPM->x.rows() == nonSymmetricIPM->c.rows());

        return nonSymmetricIPM;
    }

    IPMDouble calc_step_length_predictor() {
        IPMDouble const epsilon = .5;
        IPMDouble const eta = _large_neighborhood * pow(epsilon, _num_corrector_steps);
        IPMDouble const k_x = eta + sqrt(2 * eta * eta + _barrier->concordance_parameter(x) + 1);
        return _param_step_length_predictor / k_x;
    }

    int run_solver();

    enum Termination {
        SUCCESS,
        FAILURE
    };

    inline IPMDouble primal_error() {
        return (A * x - tau * b).norm() / _err_consts.primal;
    }

    inline IPMDouble dual_error() {
        return (A.transpose() * y + s - tau * c).norm() / _err_consts.dual;
    }

    inline IPMDouble primal_error_rescaled() {
        return (A * x / tau - b).norm();
    }

    inline IPMDouble dual_error_rescaled() {
        return (A.transpose() * y / tau + s / tau - c).norm();
    }

    inline IPMDouble duality_gap() {
        return x.dot(s) / _barrier->concordance_parameter(x);

    }

    inline IPMDouble complementarity() {
        return abs(c.dot(x) - b.dot(y)) / _err_consts.complementary;
    }

    bool verify_solution(IPMDouble precision = 10e-5) {
        if (not _barrier->in_interior(x)) {
            return false;
        }
        if (primal_error_rescaled() > precision) {
            return false;
        }
        if (dual_error_rescaled() > precision) {
            return false;
        }

        //TODO: Add duality gap check. In general we don't have strict duality as exploited below.
        if (duality_gap() > precision) {
            return false;
        }
        return true;
    }

    Solution<IPMDouble> get_solution() {
        Solution<IPMDouble> sol;
        assert(tau != 0);
        sol.x = x / tau;
        sol.s = s / tau;
        sol.centrality = centrality();
        sol.gap = mu();
        return sol;
    }

    std::shared_ptr<spdlog::logger> _logger;
    std::shared_ptr<spdlog::logger> _benchmark_logger;


    //TODO:Figure out how to make these variables while being able to compile it.

    Matrix A;
    Vector b;
    Vector c;

    Eigen::SparseMatrix<IPMDouble> A_sparse;

    //Matrix whose columns are a basis of the kernel of A. Currently unused.
    Matrix _basis_ker_A;

    Vector x;
    Vector y;
    Vector s;

    pt::ptree _config;

    unsigned _num_predictor_steps = 500;
    unsigned _num_corrector_steps;

    IPMDouble _param_step_length_predictor = 0.02;

    IPMDouble _step_length_predictor;
    IPMDouble _step_length_corrector;

    IPMDouble _epsilon = 10e-5;

    unsigned _total_num_line_steps;

    //Large neighborhood
    IPMDouble _large_neighborhood;
    //Small neighborhood
    IPMDouble _small_neighborhood;

    bool _check_centrality_in_every_segment = true;

    bool _type_cast_if_unsuccessful = true;

    bool _use_line_search = true;

    IPMDouble kappa, tau;

    Vector _last_predictor_direction;

    ErrorConstants<IPMDouble> _err_consts;

    void initialize();

private:

    cxxtimer::Timer _predictor_timer;
    cxxtimer::Timer _corrector_timer;
    cxxtimer::Timer _andersen_sys_timer;
    cxxtimer::Timer _centrality_timer;

    std::vector<cxxtimer::Timer> _custom_timers;

    cxxtimer::Timer _general_method_timer;
    cxxtimer::Timer _total_runtime_timer;


    void set_configuration_variables();

    bool terminate();

    bool terminate_successfully_wrapper();

    bool terminate_successfully();

    bool terminate_infeasible_wrapper();

    bool terminate_infeasible();


    Vector _stored_x_centrality;
    Vector _stored_s_centrality;
    IPMDouble _stored_centrality_error;


    IPMDouble mu();

    Vector psi(IPMDouble t);

    LHSCB<IPMDouble> *_barrier;

    std::vector<std::pair<Vector, Vector> >
    solve_andersen_andersen_subsystem(std::vector<std::pair<Vector, Vector> > &);

    Vector andersen_andersen_solve(Vector const);

    Vector solve(Matrix &, Vector const);

    Vector create_predictor_RHS();

    Vector create_corrector_RHS();

    Vector solve_predictor_system();

    Vector solve_corrector_system();

    IPMDouble centrality();


    void print_status();

    void apply_update(Vector concat);

    Vector build_update_vector();

    void test_gradient();

    void test_hessian();
};

#include "NonSymmetricIPM.hpp"

#endif //SOS_NONSYMMETRICIPM_H
