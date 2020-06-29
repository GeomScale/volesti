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
//    assert(r1.rows() == A.rows());
//    assert(r2.rows() == A.cols());

    Matrix mu_H_x = mu() * _barrier->hessian(x);

    //TODO: Can A(LL^\top)A^\top be computed faster?
    //TODO: Inverse maintenance, in particular for corrector steps.

    Matrix normalized_inverse_hessian = _barrier->inverse_hessian(x) / mu();
    Matrix A_H_inv_A_top = -A * normalized_inverse_hessian * A.transpose();

    std::vector<std::pair<Vector, Vector> > results;
    for (int i = 0; i < v.size(); i++) {
        auto r1 = v[i].first;
        auto r2 = v[i].second;
        auto new_s = solve(A_H_inv_A_top, A * normalized_inverse_hessian * r2 - r1);
        //TODO: might be more stable to formulate next line as linear system solve instead of using the inverse.
        auto new_t = normalized_inverse_hessian * (r2 + A.transpose() * new_s);
        results.emplace_back(std::pair<Vector, Vector>(new_s, new_t));
    }
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
    _specific_method_timer.start();
    std::vector<std::pair<Vector, Vector> > ret = solve_andersen_andersen_subsystem(new_rhs_vectors);
    _specific_method_timer.stop();

    Vector new_p = ret[0].first;
    Vector new_q = ret[0].second;

    Vector new_u = ret[1].first;
    Vector new_v = ret[1].second;

    Matrix K = Matrix::Zero(m + n, m + n);
    K.block(0, m, m, n) = A;
    K.block(m, 0, n, m) = -A.transpose();
    K.block(m, m, n, n) = mu_H_x;

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

    _logger->trace("{}", new_sol_pq.segment(m, n).transpose());
    _logger->trace("{}", new_q.transpose());

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

    _logger->info("G: \n {}", G);

    Matrix M = create_skajaa_ye_matrix();
    _logger->debug("\n M has dimension ({}, {})", M.rows(), M.cols());

    _logger->info("b: {}", b.transpose());
    _logger->info("c: {}", c.transpose());
    _logger->info("gradient: {}", _barrier->gradient(x).transpose());

    print();

    for (int pred_iteration = 0; pred_iteration < _num_predictor_steps; ++pred_iteration) {
        _logger->debug("Begin predictor iteration {}", pred_iteration);
        _predictor_timer.start();
        if (terminate()) {
            _logger->debug("terminate successfully with required proximity.");
            break;
        }
//        test_hessian();
        //TODO: rewrite code. Currently concatenation is not very smooth. Pass values by reference to concat.
        create_skajaa_ye_matrix();
        //Begin Debug content

        assert(centrality() < _beta);
        //End debug content

        _logger->trace("Begin solving predictor system");
        Vector predictor_direction = solve_predictor_system();
        _logger->trace("Finished solving predictor system");

        _logger->debug("Applied RHS for predictor direction is {}", ((_M * predictor_direction)).transpose());

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
            while ((kappa > 0 and tau > 0 and _barrier->in_interior(x)) and centrality() < _beta and not terminate()
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

            assert(terminate() or num_line_steps > 0);
        }

        _logger->info("End of predictor step {} with {} line steps:", pred_iteration, num_line_steps);
        if (_logger->level() <= spdlog::level::info) {
            print();
        }

        _predictor_timer.stop();
        _corrector_timer.start();

        for (int corr_iteration = 0; corr_iteration < _num_corrector_steps; ++corr_iteration) {
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

            assert(kappa > 0);
            assert(tau > 0);
            assert(_barrier->in_interior(x));
            assert(centrality() < _beta);
        }

        _corrector_timer.stop();
    }
}

NonSymmetricIPM::NonSymmetricIPM(Matrix &A_, Vector &b_, Vector &c_, LHSCB *barrier_) :
        A(A_), b(b_), c(c_), _barrier(barrier_), kappa(1.), tau(1.0) {
    y = Matrix::Zero(A.rows(), 1);

    _logger = spdlog::get("NonSymmetricIPM");
    if (_logger == nullptr) {
        _logger = spdlog::stdout_color_mt("NonSymmetricIPM");
        _logger->set_level(spdlog::level::debug);
    }
    initialize();
}

void NonSymmetricIPM::initialize() {
    x = _barrier->initialize_x();
    s = _barrier->initialize_s();

    std::cout << "x is initialized as: " << x.transpose() << std::endl;
    std::cout << "s is initialized as: " << s.transpose() << std::endl;

    Matrix QR = A.transpose().householderQr().householderQ();

    _basis_ker_A = QR.block(0, A.rows(), QR.rows(), QR.cols() - A.rows());

    _logger->info("Matrix A is: \n {}", A);
    _logger->info("Basis of ker A is: \n {}", _basis_ker_A);
    _logger->info("Check correctness: \n {}", A * _basis_ker_A);
    //Not exactly 0 for numerical reasons.
//    assert(A * _basis_ker_A == Matrix::Zero(A.rows(), A.cols() - A.rows()));

    _num_corrector_steps = 3;
    _beta = .99;
    _beta_small = 0.1;
    IPMDouble epsilon = 0.5;
    IPMDouble eta = _beta * pow(epsilon, _num_corrector_steps);
    IPMDouble k_x = eta + sqrt(2 * eta * eta + _barrier->concordance_parameter(x) + 1);

    _step_length_predictor = 0.020 / k_x;
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
    Vector rhs(first_set_of_equalities.rows() + s.rows() + 1);
    Vector second_set_of_equalities = -psi(mu());
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
    _logger->debug("System RHS for predictor step is {}", rhs.transpose());
    Vector andersen_direction = andersen_andersen_solve(rhs);
    return andersen_direction;
}

Vector NonSymmetricIPM::solve_corrector_system() {
    Vector rhs = create_corrector_RHS();
    _logger->debug("Psi is {}", psi(mu()).transpose());
    _logger->debug("Barrier is {}", _barrier->gradient(x).transpose());
    _logger->debug("s is {}", s.transpose());
    _logger->debug("System RHS for corrector step is {}", rhs.transpose());
    _logger->trace("Corrector matrix is: \n{}", _M);

    auto andersen_dir = andersen_andersen_solve(rhs);

    //Test in original system

    create_skajaa_ye_matrix();
    Matrix aux_sy(2, rhs.rows());
    aux_sy.block(0, 0, 1, rhs.rows()) = rhs.transpose();
    aux_sy.block(1, 0, 1, rhs.rows()) = (_M * andersen_dir).transpose();

    _logger->trace("Compare solutions: \n{}", aux_sy);

    //show that new psi vector is pretty close to 0


    auto delta_s = andersen_dir.segment(A.rows() + A.cols() + 1, A.cols());
    auto delta_x = andersen_dir.segment(A.rows(), A.cols());

//    std::cout << "Check validity of solution: "
//    << (mu() * _barrier->hessian(x) * delta_x + delta_x + psi(mu()).segment(0,delta_x.rows())).transpose() << std::endl;

//    std::cout << "Validity of big solution"
//    << (_M.block(A.cols() + A.rows() + 1, 0, A.cols(), _M.cols())
//    * andersen_dir + psi(mu()).segment(0,delta_x.rows())).transpose()
//    << std::endl;

    Vector tmp(s.rows() + 1);
    tmp.block(0, 0, s.rows(), 1) = s + mu() * _barrier->gradient(x);
    tmp(s.rows()) = kappa - mu() / tau;
//    std::cout << "Validity check of psi: " << (psi(mu()) - tmp).transpose() << std::endl;

//    std::cout << "Current vectors are x: \n" << x.transpose() << std::endl << "delta x \n" << delta_x.transpose() << std::endl;
//    std::cout << "Current vectors are s: \n" << s.transpose() << std::endl << "delta s \n" << delta_s.transpose() << std::endl;

    Vector new_psi_vector = s + delta_s + mu() * _barrier->gradient(x + delta_x);
    Matrix aux(2, new_psi_vector.rows());
    Vector calculated_psi_vector = s + delta_s
                                   + mu() * (_barrier->gradient(x) + _barrier->hessian(x) * delta_x);
    aux.block(0, 0, 1, new_psi_vector.rows()) = new_psi_vector.transpose();
    aux.block(1, 0, 1, new_psi_vector.rows()) = calculated_psi_vector.transpose();

//    Vector dir = solve(_M, rhs);
    return andersen_dir;
}

void NonSymmetricIPM::print() {

    _logger->info("alpha_pred: {}", _step_length_predictor);
    _logger->info("alpha_corrector: {}", _step_length_corrector);

    Matrix aux(2, x.rows());
    aux.block(0, 0, 1, x.rows()) = x.transpose();
    aux.block(1, 0, 1, s.rows()) = s.transpose();

    _logger->debug("Current primal/dual x, s pair :\n{}", aux);
    _logger->debug("Current primal/dual rescaled x, s pair :\n{}", aux / tau);

    _logger->info("kappa: {}, tau: {}", kappa, tau);
    _logger->info("mu: {}", mu() / (tau * tau));
    _logger->info("mu unscaled: {}", mu());
    _logger->info("centrality error: {}", centrality());
    _logger->info("duality gap: {}", kappa / tau);
    _logger->info("centrality error induced by kappa/tau: {}", (tau * kappa - mu()) / mu());
    _logger->trace("last predictor direction: {}", _last_predictor_direction.transpose());
    _logger->info("primal feasibility error: {}", (A * x / tau - b).norm());
    _logger->info("primal feasibility error unscaled: {}", (A * x - b * tau).norm());
    _logger->info("dual feasibility error: {}", (A.transpose() * y / tau + s / tau - c).norm());
    _logger->info("dual feasibility error unscaled: {}", (A.transpose() * y + s - c * tau).norm());
    _logger->info("Total elapsed predictor time: {} seconds.", _predictor_timer.count<std::chrono::seconds>());
    _logger->info("Total elapsed corrector time: {} seconds.", _corrector_timer.count<std::chrono::seconds>());
    _logger->info("Total elapsed andersen sys solve time: {} seconds.",
                  _andersen_sys_timer.count<std::chrono::seconds>());
    _logger->info("Total elapsed calc centrality time: {} seconds.", _centrality_timer.count<std::chrono::seconds>());
    _logger->debug("Total elapsed time checking if in interior: {} seconds.",_barrier->_in_interior_timer.count<std::chrono::seconds>());

    _logger->info("Total elapsed time in general method: {} seconds.",_general_method_timer.count<std::chrono::seconds>());
    _logger->info("Total elapsed time in specific method: {} seconds.",_specific_method_timer.count<std::chrono::seconds>());

    _logger->info("--------------------------------------------------------------------------------------");

}

bool NonSymmetricIPM::terminate() {
    //Duality
    if (x.dot(s) > _epsilon * tau * tau) {
        return false;
    }
    //Primal feasibility
    if ((A * x / tau - b).norm() > _epsilon) {
        return false;
    }
    //Dual feasibility
    if ((A.transpose() * y / tau + s / tau - c).norm() > _epsilon) {
        return false;
    }
    return true;
}

IPMDouble NonSymmetricIPM::centrality() {

    _centrality_timer.start();

    IPMDouble mu_d = mu();
    Vector psi_vec = psi(mu_d);

//    auto hess = _barrier->hessian(x);
    Eigen::LLT<Matrix> LLT = _barrier->llt(x);

    if (LLT.info() == Eigen::NumericalIssue) {
        std::cout << "Issue in LLT decomposition";
        assert(false);
    }

    Vector LLT_sol = LLT.matrixL().solve(psi_vec.segment(0, psi_vec.rows() - 1));
    IPMDouble tau_kappa_entry = tau * psi_vec.segment(psi_vec.rows() - 1, 1).sum();

    Vector err_L(psi_vec.rows());
    err_L << LLT_sol, tau_kappa_entry;

    IPMDouble centr_err_L = err_L.norm() / mu_d;

//    Matrix full_hessian = Matrix::Zero(A.cols() + 1, A.cols() + 1);
//    full_hessian.block(0, 0, A.cols(), A.cols()) = _barrier->hessian(x);
//    full_hessian.block(A.cols(), A.cols(), 1, 1) = 1 / (tau * tau) * Matrix::Identity(1, 1);
//    auto sys_solve = solve(full_hessian, psi_vec);
//    Double centr_err2 = boost::multiprecision::sqrt(boost::multiprecision::abs((psi_vec.transpose() * sys_solve).sum())) / mu_d;

//    assert(centr_err2 - centr_err < 10e-5);

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

