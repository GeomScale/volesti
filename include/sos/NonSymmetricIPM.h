//
// Created by Bento Natura on 07/06/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_NONSYMMETRICIPM_H
#define NONSYMMETRICCONICOPTIMIZATION_NONSYMMETRICIPM_H

#include "LHSCB.h"
#include <iostream>
#include <cxxtimer.hpp>


class Instance {
public:
    Constraints constraints;
    LHSCB *barrier;
};

class DirectionDecomposition {
public:
    DirectionDecomposition(Vector v, unsigned const n, unsigned const m) {
        assert(v.rows() == m + n + 1 + n + 1);
        y = v.block(0, 0, m, 1);
        x = v.block(m, 0, n, 1);
        tau = v.block(m + n, 0, 1, 1).sum();
        s = v.block(m + n + 1, 0, n, 1);
        kappa = v.block(m + n + 1 + n, 0, 1, 1).sum();
    }

    friend std::ostream &operator<<(std::ostream &os, const DirectionDecomposition &dir);

    Vector y, x, s;
    IPMDouble kappa, tau;
};

class NonSymmetricIPM {

public:

    NonSymmetricIPM(Matrix &A_, Vector &b_, Vector &c_, LHSCB *barrier_);

    NonSymmetricIPM(Instance &instance) : NonSymmetricIPM(instance.constraints.A,
                                                          instance.constraints.b, instance.constraints.c,
                                                          instance.barrier) {
        _logger->set_level(spdlog::level::info);
    };

    void run_solver();

    inline IPMDouble primal_error() {
        return (A * x / tau - b).norm();
    }

    inline IPMDouble dual_error() {
        return (A.transpose() * y / tau + s / tau - c).norm();
    }

    inline IPMDouble duality_gap(){
        return x.dot(s) / _barrier->concordance_parameter(x);
    }

    bool verify_solution(IPMDouble precision = 10e-5) {
        if (not _barrier->in_interior(x)) {
            return false;
        }
        if (primal_error() > precision) {
            return false;
        }
        if (dual_error() > precision) {
            return false;
        }

        //TODO: Add duality gap check. In general we don't have strict duality as exploited below.
        if (duality_gap() > precision) {
            return false;
        }

        return true;
    }

    Solution get_solution() {
        Solution sol;
        assert(tau != 0);
        sol.x = x / tau;
        sol.s = s / tau;
        sol.centrality = centrality();
        sol.gap = mu();
        return sol;
    }

    std::shared_ptr<spdlog::logger> _logger;

private:

    Matrix A;
    Vector b;
    Vector c;

    Matrix _basis_ker_A;

    Vector x;
    Vector y;
    Vector s;

    Matrix _M;
    Matrix _G;

    unsigned _num_predictor_steps = 500;
    unsigned _num_corrector_steps;
    IPMDouble _step_length_predictor;
    IPMDouble _step_length_corrector;

    IPMDouble _epsilon = 10e-5;

    cxxtimer::Timer _predictor_timer;
    cxxtimer::Timer _corrector_timer;
    cxxtimer::Timer _andersen_sys_timer;
    cxxtimer::Timer _centrality_timer;

    cxxtimer::Timer _general_method_timer;
    cxxtimer::Timer _specific_method_timer;

    bool terminate();

    bool _use_line_search = true;

    IPMDouble kappa, tau;

    //Large neighborhood
    IPMDouble _beta;
    //Small neighborhood
    IPMDouble _beta_small;

    IPMDouble mu();

    Vector psi(IPMDouble t);

    LHSCB *_barrier;

    Matrix create_skajaa_ye_matrix();

    std::vector<std::pair<Vector, Vector> > solve_andersen_andersen_subsystem(std::vector<std::pair<Vector, Vector> > &);

    Vector andersen_andersen_solve(Vector const);

    Vector solve(Matrix &, Vector const);

    Matrix create_matrix_G();

    Vector create_predictor_RHS();

    Vector create_corrector_RHS();

    Vector solve_predictor_system();

    Vector solve_corrector_system();

    IPMDouble centrality();


    void print();

    void apply_update(Vector concat);

    Vector build_update_vector();

    Vector _last_predictor_direction;

    void initialize();

    void test_gradient();

    void test_hessian();
};


#endif //NONSYMMETRICCONICOPTIMIZATION_NONSYMMETRICIPM_H
