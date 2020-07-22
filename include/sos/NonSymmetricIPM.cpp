// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE filee

#include <spdlog/sinks/basic_file_sink.h>
#include "NonSymmetricIPM.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/fmt/ostr.h"


std::ostream &operator<<(std::ostream &os, const DirectionDecomposition &dir) {
    os << "x " << dir.x.transpose() << std::endl;
    os << "s: " << dir.s.transpose() << std::endl;
    os << "y: " << dir.y.transpose() << std::endl;
    os << "kappa: " << dir.kappa << ", tau: " << dir.tau << std::endl;
    return os;
}

Vector NonSymmetricIPM::solve(Matrix &M_, Vector const v_) {
    Vector sol = M_.colPivHouseholderQr().solve(v_);
    return sol;
}

Matrix NonSymmetricIPM::create_matrix_G() {
    const int m = A.rows();
    const int n = A.cols();

    assert(b.rows() == m);
    assert(c.rows() == n);

    const int dim = m + n + 1;
    Matrix G = Matrix::Zero(dim, dim);

    G.block(0, m, m, n) = A;
    G.block(0, m + n, m, 1) = -b;
    G.block(m, 0, n, m) = -A.transpose();
    G.block(m, m + n, n, 1) = c;
    G.block(m + n, 0, 1, m) = b.transpose();
    G.block(m + n, m, 1, n) = -c.transpose();

    _G = G;
    return G;
}

//TODO:For sparse systems we might create a dense matrix in the LLS decomposition. Implement separate solver for this case

std::vector<std::pair<Vector, Vector> > NonSymmetricIPM::solve_andersen_andersen_subsystem(
        std::vector<std::pair<Vector, Vector> > &v) {


    _test_timers[3].start();
    Eigen::LLT<Matrix> LLT = _barrier->llt(x);
    LLT.matrixL().toDenseMatrix().eval();
    _test_timers[3].stop();

    //TODO: Can A(LL^\top)A^\top be computed faster?
    //TODO: Inverse maintenance, in particular for corrector steps.

    _test_timers[8].start();
    _test_timers[4].start();
//    Matrix normalized_inverse_hessian = _barrier->inverse_hessian(x) / mu();
//    normalized_inverse_hessian.eval();

    //TODO: Turn into blockwise operation
//    Matrix A_H_inv = A * normalized_inverse_hessian;
    //TODO: double transposition. Figure out how to multiply solve from RHS.
//    Matrix A_H_inv = LLT.solve(A.transpose()).transpose() / mu();
    Matrix A_H_inv = _barrier->llt_solve(x, A.transpose()).transpose() / mu();
//    A_H_inv.eval();
    A_H_inv.eval();
    _test_timers[4].stop();
    //Note: Below is another valid version to compute A_H_inv_alt, but it appears to be slower.
//    Matrix L_inv = LLT.matrixL().solve(Matrix::Identity(x.rows(),x.rows()));
//    Matrix A_H_inv_alt = A_sparse * L_inv.transpose() *  L_inv.triangularView<Eigen::Lower>() / mu();
//    Matrix A_H_inv_alt = LLT.matrixL().solve((A_sparse * L_inv.transpose()).transpose()).transpose() / mu();
//    A_H_inv_alt.eval();

//    std::cout << "A_H_inv is \n" << A_H_inv << std::endl;
//    std::cout << "A_H_inv_alt is \n" << A_H_inv_alt << std::endl;

    _test_timers[5].start();

    //TODO: Use better method to sparsify A.
    Matrix A_H_inv_A_top = A_H_inv * A_sparse.transpose();
    A_H_inv_A_top.eval();

    _test_timers[5].stop();
    _test_timers[6].start();
    Eigen::LLT<Matrix> A_H_inv_A_top_LLT = A_H_inv_A_top.llt();
    A_H_inv_A_top.llt().matrixL().toDenseMatrix().eval();
    _test_timers[6].stop();
    _test_timers[8].stop();
    _test_timers[9].start();

    std::vector<std::pair<Vector, Vector> > results;
    //TODO: might be possible to solve these system in batches! (But currently runtime of the loop is clearly dominated by previous parts)
    //TODO: check whether Conjugate Gradient Method solves Normal Equations more efficiently.
    for (unsigned i = 0; i < v.size(); i++) {
        Vector &r1 = v[i].first;
        Vector &r2 = v[i].second;
        Vector new_s_intermediate = A_H_inv_A_top_LLT.matrixL().solve(-(A_H_inv * r2 - r1));
        Vector new_s = A_H_inv_A_top_LLT.matrixU().solve(new_s_intermediate);
        //TODO: might be more stable to formulate next line as linear system solve instead of using the inverse.
        Vector new_t_intermediate = LLT.matrixL().solve((r2 + A.transpose() * new_s) / mu());
        Vector new_t = LLT.matrixU().solve(new_t_intermediate);
        results.emplace_back(std::pair<Vector, Vector>(new_s, new_t));
    }
    _test_timers[9].stop();
    return results;
}

Vector NonSymmetricIPM::andersen_andersen_solve(Vector const rhs) {

    _andersen_sys_timer.start();

    const int m = A.rows();
    const int n = A.cols();

    //TODO: figure out whether rescaling makes sense.
    Vector const r_p = rhs.segment(0, m);
    Vector const r_d = rhs.segment(m, n);
    Vector const r_g = rhs.segment(m + n, 1);
    Vector const r_xs = rhs.segment(m + n + 1, n);
    Vector const r_tk = rhs.segment(m + n + 1 + n, 1);

    Matrix mu_H_x = mu() * _barrier->hessian(x);
    IPMDouble mu_H_tau = mu() / (tau * tau);

    std::vector<std::pair<Vector, Vector> > new_rhs_vectors;
    new_rhs_vectors.emplace_back(std::pair<Vector, Vector>(b, -c));
    new_rhs_vectors.emplace_back(std::pair<Vector, Vector>(r_p, r_d + r_xs));

    //TODO: runtime test of both methods to solve the subsystem. Currently both are performed, but eventually
    // we'll stick with the faster one.
    std::vector<std::pair<Vector, Vector> > ret = solve_andersen_andersen_subsystem(new_rhs_vectors);

    Vector new_p = ret[0].first;
    Vector new_q = ret[0].second;

    Vector new_u = ret[1].first;
    Vector new_v = ret[1].second;

//    Matrix K = Matrix::Zero(m + n, m + n);
//    K.block(0, m, m, n) = A;
//    K.block(m, 0, n, m) = -A.transpose();
//    K.block(m, m, n, n) = mu_H_x;

    Vector rhs_pq(m + n);
    rhs_pq << b, -c;

    Vector aux_new_pq(new_p.rows() + new_q.rows());
    aux_new_pq << new_p, new_q;

//    _logger->info("Orig: {}, \n New {}", pq.segment(0, m).transpose(),
//                   aux_new_pq.transpose());
//    _logger->info("Diff: {}", (aux_new_pq - pq).norm());

    Vector new_sol_pq(m + n);
    new_sol_pq << new_p, new_q;

    Vector new_sol_uv(m + n);
    new_sol_uv << new_u, new_v;

    SPDLOG_TRACE("{}", new_sol_pq.segment(m, n).transpose());
    SPDLOG_TRACE("{}", new_q.transpose());

    Vector rhs_uv(m + n);
    rhs_uv << r_p, r_d + r_xs;

    IPMDouble d_tau = (r_g.sum() + r_tk.sum() - rhs_pq.dot(new_sol_uv)) / (mu() / (tau * tau) + rhs_pq.dot(new_sol_pq));

    Vector d_yx = new_sol_uv + d_tau * new_sol_pq;
    Vector d_x = d_yx.segment(m, n);
    Vector d_s = r_xs - mu_H_x * d_x;
    IPMDouble d_kappa = r_tk.sum() - mu_H_tau * d_tau;

    Vector d_tau_vec(1);
    d_tau_vec(0) = d_tau;
    Vector d_kappa_vec(1);
    d_kappa_vec(0) = d_kappa;

    Vector sol(m + n + 1 + n + 1);
    sol << d_yx, d_tau_vec, d_s, d_kappa_vec;

    _andersen_sys_timer.stop();

    //Check for dual error in this solution:

    Vector dual_err = -A.transpose() * d_yx.segment(0, y.rows()) + d_tau * c - d_s - r_d;
    IPMDouble rel_dual_error = dual_err.norm() / r_d.norm();

    _logger->debug("relative dual error is {}", rel_dual_error);
    if (rel_dual_error > 10e-2) {
        _logger->warn("relative dual error is {}", rel_dual_error);
    }

    return sol;
}

Matrix NonSymmetricIPM::create_skajaa_ye_matrix() {
    _logger->trace("Begin creating skajaa ye matrix");
    const int m = A.rows();
    const int n = A.cols();
    Matrix G = create_matrix_G();
    const int dim = G.rows() + n + 1;
    Matrix M = Matrix::Zero(dim, dim);
    M.block(0, 0, G.rows(), G.cols()) = G;
    M.block(m, G.cols(), n + 1, n + 1) = -Matrix::Identity(n + 1, n + 1);

    //second set of inequalities
    M.block(G.rows(), m, n, n) = mu() * _barrier->hessian(x);
    M.block(G.rows(), m + n + 1, n, n) = Matrix::Identity(n, n);

    //last row entry for kappa
    M.block(G.rows() + n, m + n + 1 + n, 1, 1) = Matrix::Identity(1, 1);

    //last row entry for tau
    M.block(G.rows() + n, m + n, 1, 1) =
            mu() / (tau * tau) * Matrix::Identity(1, 1);

    _M = M;

    _logger->trace("Hessian for vector is: \n{}", _barrier->hessian(x));
    _logger->trace("Constraint matrix for system: \n{}", M);

    _logger->trace("Finished creating skajaa ye matrix");
    return M;
}


Vector NonSymmetricIPM::build_update_vector() {
    Vector concat(y.rows() + x.rows() + 1 + s.rows() + 1);
    concat << y, x, tau * Matrix::Identity(1, 1), s, kappa * Matrix::Identity(1, 1);
    return concat;
}

void NonSymmetricIPM::apply_update(Vector concat) {
    y = concat.block(0, 0, y.rows(), 1);
    x = concat.block(y.rows(), 0, x.rows(), 1);
    tau = concat.block(y.rows() + x.rows(), 0, 1, 1).sum();
    s = concat.block(y.rows() + x.rows() + 1, 0, s.rows(), 1);
    kappa = concat.block(y.rows() + x.rows() + 1 + s.rows(), 0, 1, 1).sum();
}

void NonSymmetricIPM::run_solver() {

    Matrix G = create_matrix_G();

    _logger->debug("G: \n {}", G);

    Matrix M = create_skajaa_ye_matrix();
    _logger->debug("\n M has dimension ({}, {})", M.rows(), M.cols());

    _logger->debug("b: {}", b.transpose());
    _logger->debug("c: {}", c.transpose());
    _logger->debug("gradient: {}", _barrier->gradient(x).transpose());

    print();

    _total_num_line_steps = 0;

    _total_runtime_timer.start();
    for (unsigned pred_iteration = 0; pred_iteration < _num_predictor_steps; ++pred_iteration) {
        _logger->debug("Begin predictor iteration {}", pred_iteration);
        _predictor_timer.start();
        if (terminate_successfully_wrapper()) {
            _logger->info("Interior point method terminated successfully with required proximity.");
            break;
        }

        if (terminate_infeasible_wrapper()) {
            _logger->info("Interior point method terminated with infeasible solution.");
            break;
        }
//        test_hessian();
        //TODO: rewrite code. Currently concatenation is not very smooth. Pass values by reference to concat.
//        create_skajaa_ye_matrix();

        assert(centrality() < _beta);

        _logger->trace("Begin solving predictor system");
        Vector predictor_direction = solve_predictor_system();
        _logger->trace("Finished solving predictor system");

//        _logger->debug("Applied RHS for predictor direction is {}", ((_M * predictor_direction)).transpose());

        DirectionDecomposition dir(_step_length_predictor * predictor_direction, x.rows(), y.rows());
        _logger->debug("Predictor direction is \n {}", dir);
        Vector concat = build_update_vector();
        unsigned num_line_steps = 0;
        if (not _use_line_search) {
            concat += _step_length_predictor * predictor_direction;
            apply_update(concat);

        } else {
            //find longest step length such that we remain in beta environment.
            //TODO: find sophisticated way of computing this efficiently. Currently we do simple repeated squaring,
            // irrespective of the barrier function.

            concat += _step_length_predictor * predictor_direction;
            apply_update(concat);
            Vector fallback_vec = concat;
            while (kappa > 0 and tau > 0 and _barrier->in_interior(x) and centrality() < _beta and not terminate()
                   and pow(2, num_line_steps) * _step_length_predictor <= .5) {
                _logger->debug("Another iteration in line search...");
//                print();
                fallback_vec = concat;
                concat += pow(2, num_line_steps) * _step_length_predictor * predictor_direction;
                num_line_steps++;
                apply_update(concat);
            }

            _logger->debug("Reason for stopping: ");
            if (not _barrier->in_interior(x)) {
                _logger->debug("Not in interior");
            } else if (centrality() >= _beta) {
                _logger->debug("Centrality {} is worse than neighborhood {}", centrality(), _beta);
            }
            _logger->debug("Applied {} line steps in iteration {}", num_line_steps, pred_iteration);
            apply_update(fallback_vec);

            //In this case we perform more corrector steps instead.
//            assert(terminate_successfully() or num_line_steps > 0);
        }

        _total_num_line_steps += num_line_steps;
        _logger->info("End of predictor step {} with {} line steps and total num line steps {}:",
                      pred_iteration, num_line_steps, _total_num_line_steps);
        if (_logger->level() <= spdlog::level::info) {
            print();
        }

        _predictor_timer.stop();
        _corrector_timer.start();

        for (unsigned corr_iteration = 0; corr_iteration < _num_corrector_steps; ++corr_iteration) {
            //TODO: figure out if this is already as expensive as running another corrector step. (probably not, as the crrent value can be used for next predictor step if true).
            if (centrality() < _beta_small) {
                break;
            }
            create_skajaa_ye_matrix();

            Vector psi_full = psi(mu());

            //TODO: Use Woodburry matrix identity for corrector step.

            Vector corrector_direction = solve_corrector_system();

            _logger->debug("RHS for corrector direction is {}", ((_M * corrector_direction)).transpose());
            DirectionDecomposition dir(corrector_direction, x.rows(), y.rows());
            _logger->debug("Corrector direction is: \n{}", dir);
            Vector concat = build_update_vector();
            concat += _step_length_corrector * corrector_direction;
            apply_update(concat);

            _logger->info("End of corrector step {} :", corr_iteration);
            if (_logger->level() <= spdlog::level::info) {
                print();
            }

            if (_logger->level() <= spdlog::level::debug) {
                assert(kappa > 0);
                assert(tau > 0);
                assert(_barrier->in_interior(x));
                assert(centrality() < _beta);
            }
        }

        _corrector_timer.stop();
    }
    _total_runtime_timer.stop();
}

NonSymmetricIPM::NonSymmetricIPM(Matrix &A_, Vector &b_, Vector &c_, LHSCB *barrier_) :
        A(A_), b(b_), c(c_), kappa(1.), tau(1.), _barrier(barrier_) {
    y = Matrix::Zero(A.rows(), 1);

    //TODO: use proper tolerance / reference.
    A_sparse = A.sparseView(IPMDouble(1e-10), 1e-10);

    _stored_x_centrality.resize(c.rows());
    _stored_s_centrality.resize(c.rows());

    _test_timers.resize(10);

    _logger = spdlog::get("NonSymmetricIPM");
    if (_logger == nullptr) {
        std::vector<spdlog::sink_ptr> sinks;
        sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
        sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/logfile.txt"));
        _logger = std::make_shared<spdlog::logger>("NonSymmetricIPM", begin(sinks), end(sinks));
        _logger->set_level(spdlog::level::info);

    }
    initialize();
}

void NonSymmetricIPM::initialize() {

    //Rescale instance for stability/conditioning

    IPMDouble scaling_delta_primal = 0;
    for (int i = 0; i < A.rows(); i++) {
        IPMDouble const row_ratio = (1. + abs(b(i))) / (1. + abs(A.row(i).sum()));
        scaling_delta_primal = std::max(scaling_delta_primal, row_ratio);
    }

    x = _barrier->initialize_x();
    s = _barrier->initialize_s();

    IPMDouble scaling_delta_dual = 0;
    for (int i = 0; i < s.rows(); i++) {
        IPMDouble const entry_ratio = (1 + abs(s(i))) / (1 + abs(c(i)));
        scaling_delta_dual = std::max(scaling_delta_dual, entry_ratio);
    }

    IPMDouble const scaling_delta = sqrt(scaling_delta_dual * scaling_delta_primal);

    _logger->info("Norm of x is {} and norm of s is {} before rescaling.", x.norm(), s.norm());

    x = _barrier->initialize_x(scaling_delta);
    s = _barrier->initialize_s(scaling_delta);

    _logger->debug("Rescaled initial point by {}", scaling_delta);
    _logger->debug("Norm of c is {}, Norm of b is {} and norm of A is {} ", c.norm(), b.norm(), A.norm());

    assert(x.rows() == _barrier->getNumVariables());

    _logger->debug("Norm of x is initialized as: {}", x.norm());
    _logger->debug("Norm of s is initialized as: {}", s.norm());

    //TODO: turn into sparse operation based on n number of nonzeros.
    Matrix QR = A.transpose().householderQr().householderQ();

    Eigen::SparseQR<Eigen::SparseMatrix<IPMDouble>, Eigen::COLAMDOrdering<int> > QR_sparse;

    Eigen::SparseMatrix<IPMDouble> tmp_sparse = A_sparse.transpose();
    tmp_sparse.makeCompressed();
    QR_sparse.compute(tmp_sparse);


//    _basis_ker_A = QR.block(0, A.rows(), QR.rows(), QR.cols() - A.rows());
    _basis_ker_A = Matrix(QR_sparse.matrixQ()).block(0, A.rows(), QR_sparse.rows(), QR_sparse.cols() - A.rows());

    _logger->trace("Matrix A is: \n {}", A);
    _logger->trace("Basis of ker A is: \n {}", _basis_ker_A);
    _logger->trace("Check correctness: \n {}", A * _basis_ker_A);
    //Not exactly 0 for numerical reasons.
//    assert(A * _basis_ker_A == Matrix::Zero(A.rows(), A.cols() - A.rows()));

    _num_corrector_steps = 3;
    _beta = .99;
    _beta_small = 0.1;

    _step_length_predictor = calc_step_length_predictor();
    //TODO: Figure out what good steplength is. Even .5 step-length led to issues.
    _step_length_corrector = 1;
}

Vector NonSymmetricIPM::create_predictor_RHS() {
    Vector v(y.rows() + x.rows() + 1 + s.rows() + 1);
    v << y, x, tau * Matrix::Identity(1, 1), s, kappa * Matrix::Identity(1, 1);
    DirectionDecomposition cur_sol(v, x.rows(), y.rows());
    _logger->debug("Current vector is: \n{}", cur_sol);
    const int sub_rows = _G.rows();
    const int sub_cols = _M.cols();
    Vector first_set_of_equalities = -_M.block(0, 0, sub_rows, sub_cols) * v;
    _logger->debug("First set of equalities RHS for predictor: {} ", first_set_of_equalities.transpose());
    Vector rhs(first_set_of_equalities.rows() + s.rows() + 1);
    rhs << first_set_of_equalities, -s, -kappa * Matrix::Identity(1, 1);
    return rhs;
}

Vector NonSymmetricIPM::create_corrector_RHS() {
    Vector first_set_of_equalities = Matrix::Zero(_G.rows(), 1);
    Vector second_set_of_equalities = -psi(mu());
    Vector rhs(first_set_of_equalities.rows() + s.rows() + 1);
    rhs << first_set_of_equalities, second_set_of_equalities;
    _logger->debug("Corrector RHS is \n {}", rhs.transpose());
    return rhs;
};

IPMDouble NonSymmetricIPM::mu() {
    return (x.dot(s) + (tau * kappa)) / (_barrier->concordance_parameter(x) + 1);
};

Vector NonSymmetricIPM::psi(IPMDouble t) {
    Vector v(s.rows());
    v = s + t * _barrier->gradient(x);
    Vector v_aux = (kappa - t / tau) * Matrix::Identity(1, 1);
    Vector res(s.rows() + 1);
    res << v, v_aux;
    return res;
}

Vector NonSymmetricIPM::solve_predictor_system() {
    Vector rhs = create_predictor_RHS();
    SPDLOG_LOGGER_DEBUG(_logger, "System RHS for predictor step is {}", rhs.transpose());
    Vector andersen_direction = andersen_andersen_solve(rhs);
    return andersen_direction;
}

Vector NonSymmetricIPM::solve_corrector_system() {
    Vector rhs = create_corrector_RHS();
    SPDLOG_LOGGER_DEBUG(_logger, "Psi is {}", psi(mu()).transpose());
    SPDLOG_LOGGER_DEBUG(_logger, "Barrier is {}", _barrier->gradient(x).transpose());
    SPDLOG_LOGGER_DEBUG(_logger, "s is {}", s.transpose());
    SPDLOG_LOGGER_DEBUG(_logger, "System RHS for corrector step is {}", rhs.transpose());
    SPDLOG_LOGGER_TRACE(_logger, "Corrector matrix is: \n{}", _M);

    auto andersen_dir = andersen_andersen_solve(rhs);

    if (_logger->level() <= spdlog::level::debug) {
        //Perform iterative refinement. Should not be too expensive as costs are dominated by other stuff.
        Matrix aux_sy(2, rhs.rows());
        aux_sy.block(0, 0, 1, rhs.rows()) = rhs.transpose();
        aux_sy.block(1, 0, 1, rhs.rows()) = (_M * andersen_dir).transpose();
        _logger->debug("Compare solutions: \n{}", aux_sy);
    }
    return andersen_dir;
}

void NonSymmetricIPM::print() {

    std::string format_ = "{:<25}:{:20.2}";

    Double dummy_D;

    Double const step_length_predictor = InterpolantDoubletoIPMDouble(_step_length_predictor, dummy_D);
    Double const step_length_corrector = InterpolantDoubletoIPMDouble(_step_length_corrector, dummy_D);

    _logger->debug(format_, "alpha_predictor", step_length_predictor);
    _logger->debug(format_, "alpha_corrector", step_length_corrector);

    Matrix aux(2, x.rows());
    aux.block(0, 0, 1, x.rows()) = x.transpose();
    aux.block(1, 0, 1, s.rows()) = s.transpose();


    _logger->debug("Current primal/dual x, s pair :\n{}", aux);
    _logger->debug("Current primal/dual rescaled x, s pair :\n{}", aux / tau);

    _logger->info(format_, "kappa", InterpolantDoubletoIPMDouble(kappa, dummy_D));
    _logger->info(format_, "tau", InterpolantDoubletoIPMDouble(tau, dummy_D));

    IPMDouble mu_ipm_scaled = mu() / (tau * tau);
    Double const mu_scaled = InterpolantDoubletoIPMDouble(mu_ipm_scaled, dummy_D);
    _logger->debug(format_, "mu scaled", mu_scaled);

    IPMDouble mu_ipm = mu();
    Double const mu_ = InterpolantDoubletoIPMDouble(mu_ipm, dummy_D);
    _logger->info(format_, "mu", boost::numeric_cast<double>(mu()));

    if (_logger->level() <= spdlog::level::debug) {
        IPMDouble centrality_ipm_ = centrality();
        Double const centrality_ = InterpolantDoubletoIPMDouble(centrality_ipm_, dummy_D);
        _logger->debug(format_, "centrality error", centrality_);
    }

    IPMDouble duality_gap_ipm_ = kappa / tau;
    Double duality_gap_ = InterpolantDoubletoIPMDouble(duality_gap_ipm_, dummy_D);
    _logger->info(format_, "duality gap", duality_gap_);

    IPMDouble primal_inf_ipm_ = primal_error();
    Double primal_inf_ = InterpolantDoubletoIPMDouble(primal_inf_ipm_, dummy_D);
    _logger->info(format_, "primal infeas. ", primal_inf_);

    IPMDouble primal_inf_unscaled_ipm_ = primal_error_rescaled();
    Double primal_inf_unscaled_ = InterpolantDoubletoIPMDouble(primal_inf_unscaled_ipm_, dummy_D);
    _logger->debug(format_, "primal infeas. unscaled", primal_inf_unscaled_);

    IPMDouble dual_inf_ipm_ = dual_error();
    Double dual_inf_ = InterpolantDoubletoIPMDouble(dual_inf_ipm_, dummy_D);
    _logger->info(format_, "dual infeas.", dual_inf_);

    IPMDouble dual_inf_unscaled_ipm_ = dual_error_rescaled();
    Double dual_inf_unscaled_ = InterpolantDoubletoIPMDouble(dual_inf_unscaled_ipm_, dummy_D);
    _logger->debug(format_, "dual infeas. unscaled", dual_inf_unscaled_);

    _logger->trace("last predictor direction: {}", _last_predictor_direction.transpose());

    _logger->info(format_, "Predictor time (s)", _predictor_timer.count<std::chrono::milliseconds>() / 1000.);
    _logger->info(format_, "Corrector time (s)", _corrector_timer.count<std::chrono::milliseconds>() / 1000.);

    _logger->info(format_, "Total andersen time (s)",
                  _andersen_sys_timer.count<std::chrono::milliseconds>() / 1000.);
    _logger->info(format_, "Total runtime (s)",
                  _total_runtime_timer.count<std::chrono::milliseconds>() / 1000.);
    _logger->info(format_, "Calc centrality time (s)",
                  _centrality_timer.count<std::chrono::milliseconds>() / 1000.);
    _logger->info(format_, "Time checking interior(s)",
                  _barrier->_in_interior_timer.count<std::chrono::milliseconds>() / 1000.);


    _logger->info(format_, "Time per step: ",
                  _total_runtime_timer.count<std::chrono::milliseconds>() / (1000. * _total_num_line_steps));
    for (unsigned idx = 0; idx < _test_timers.size(); idx++) {
        std::string s = "Test timer " + std::to_string(idx);
        _logger->info(format_, s,
                      _test_timers[idx].count<std::chrono::milliseconds>() / 1000.);
    }

    _logger->info("--------------------------------------------------------------------------------------");

}

//According to Skajaa-Ye

bool NonSymmetricIPM::terminate_successfully_wrapper() {

//    Eigen::internal::set_is_malloc_allowed(true);
    bool result = terminate_successfully();
//    Eigen::internal::set_is_malloc_allowed(false);
    return result;

}

bool NonSymmetricIPM::terminate_successfully() {
    //TODO: use same criteria as in infeasilbe (i.e. scaling invariant and what is used in Skajaa-Ye)
    //Duality
    if (x.dot(s) > _epsilon * tau * tau) {
        return false;
    }
    //Primal feasibility
    if (primal_error_rescaled() > _epsilon) {
        return false;
    }

    if (primal_error() > _epsilon) {
        return false;
    }

    //Dual feasibility
    if (dual_error_rescaled() > _epsilon) {
        return false;
    }

    if (dual_error() > _epsilon) {
        return false;
    }

    return true;
}

// Termination criteria taken from Skajaa - Ye "A Homogeneous Interior-Point Algorithm for
// Nonsymmetric Convex Conic Optimization" https://web.stanford.edu/~yyye/nonsymmhsdimp.pdf page 15.

bool NonSymmetricIPM::terminate_infeasible_wrapper() {

//    Eigen::internal::set_is_malloc_allowed(true);
    bool result = terminate_infeasible();
//    Eigen::internal::set_is_malloc_allowed(false);
    return result;
}

bool NonSymmetricIPM::terminate_infeasible() {

    //TODO: Figure out if initialization scaling (delta) should influence the termination criteria.

    //Primal feasibility
    IPMDouble const ipm_1 = IPMDouble(1.);
    if ((A * x - b * tau).lpNorm<Eigen::Infinity>() >
        _epsilon * std::max(ipm_1, static_cast<IPMDouble>(A.lpNorm<Eigen::Infinity>() + b.lpNorm<Eigen::Infinity>()))) {
        return false;
    }

    //Dual feasibility
    if ((A.transpose() * y + s - c * tau).lpNorm<Eigen::Infinity>() >
        _epsilon * std::max(ipm_1, static_cast<IPMDouble>(A.lpNorm<Eigen::Infinity>() + c.lpNorm<Eigen::Infinity>()))) {
        return false;
    }

    //Duality gap
    if (abs(-c.dot(x) + b.dot(y) - kappa)
        >
        _epsilon * std::max(ipm_1, static_cast<IPMDouble>(c.lpNorm<Eigen::Infinity>() + b.lpNorm<Eigen::Infinity>()))) {
        return false;
    };

    //tiny tau
    if (tau > _epsilon * 10e-2 * std::max(ipm_1, static_cast<IPMDouble>(kappa))) {
        return false;
    }
    return true;

}

bool NonSymmetricIPM::terminate() {
    return terminate_successfully_wrapper() or terminate_infeasible_wrapper();
}


IPMDouble NonSymmetricIPM::centrality() {

    if (_stored_x_centrality == x and _stored_s_centrality == s) {
        return _stored_centrality_error;
    }
    _centrality_timer.start();

    IPMDouble mu_d = mu();

    _test_timers[0].start();
    Vector psi_vec = psi(mu_d);
    _test_timers[0].stop();

    _test_timers[2].start();
    IPMDouble tau_kappa_entry = tau * psi_vec.segment(psi_vec.rows() - 1, 1).sum();

    Vector LLT_sol = _barrier->llt_L_solve(x, psi_vec.segment(0, psi_vec.rows() - 1));
    Vector err_L(psi_vec.rows());
    err_L << LLT_sol, tau_kappa_entry;
    _test_timers[2].stop();

    IPMDouble centr_err_L = err_L.norm() / mu_d;

    _stored_x_centrality = x;
    _stored_s_centrality = s;
    _stored_centrality_error = centr_err_L;
    _centrality_timer.stop();

    return centr_err_L;
}

void NonSymmetricIPM::test_gradient() {
    //Currently the original function is not implemented.
}

void NonSymmetricIPM::test_hessian() {
    for (int i = 0; i < x.rows(); i++) {
        Vector dir = Vector::Zero(x.rows());
        dir(i) = 10e-4;
        _logger->debug("Test hessian at: {}", x.transpose());
        _logger->debug("with update {}", dir.transpose());
        Matrix M(2, x.rows());
        M.block(0, 0, 1, x.rows()) = _barrier->gradient(x + dir).transpose()
                                     - _barrier->gradient(x).transpose();
        M.block(1, 0, 1, x.rows()) = (_barrier->hessian(x) * dir).transpose();
        _logger->debug("Comparison of gradient and hessian for {}: \n{}", dir.transpose(), M);
    }
    return;
};

