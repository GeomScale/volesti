// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "../../../test/doctest.h"
#include "../utils.h"
#include "../EnvelopeProblemSDP.h"
#include "../NonSymmetricIPM.h"
#include "../barriers/LPStandardBarrier.h"
#include "../barriers/SDPStandardBarrier.h"

template<typename IPMDouble>
Constraints<IPMDouble> convert_LP_to_SDP(Matrix<IPMDouble> &A, Vector<IPMDouble> &b, Vector<IPMDouble> &c) {
    const int m = A.rows();
    const int n = A.cols();

    Constraints<IPMDouble> SDP_constraints;
    SDP_constraints.A = Matrix<IPMDouble>::Zero(m + n * n - n, n * n);

    std::cout << A << std::endl;
    for (int i = 0; i < m; ++i) {
        Matrix<IPMDouble> diag_row = A.row(i).asDiagonal();
        Eigen::Map<Matrix<IPMDouble> > diag_row_stacked(diag_row.data(), 1, n * n);
        SDP_constraints.A.block(i, 0, 1, n * n) = diag_row_stacked;
    }

    unsigned row_idx = m;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                SDP_constraints.A(row_idx, n * i + j) = 1;
                row_idx++;
            }
        }
    }

    SDP_constraints.b = Vector<IPMDouble>::Zero(SDP_constraints.A.rows());
    SDP_constraints.b.block(0, 0, m, 1) = b;

    SDP_constraints.c = Vector<IPMDouble>::Zero(SDP_constraints.A.cols());
    for (int j = 0; j < n; ++j) {
        SDP_constraints.c(j + j * n) = c(j);
    }

    return SDP_constraints;
}

template<typename IPMDouble>
bool test_sdp_solver_random_lp_formulation(const int m, const int n) {
    Matrix<IPMDouble> A = Matrix<IPMDouble>::Random(m, n);
    Vector<IPMDouble> x0 = Vector<IPMDouble>::Random(n);
    Vector<IPMDouble> c = Vector<IPMDouble>::Random(n);
    //for loop to create a feasible and bounded instance
    for (int i = 0; i < n; i++) {
        x0(i) += 1;
        c(i) += 1;
    }
    Vector<IPMDouble> b = A * x0;

    Constraints<IPMDouble> constraints = convert_LP_to_SDP<IPMDouble>(A, b, c);

    SDPStandardBarrier<IPMDouble> sdp_barrier(n);

    NonSymmetricIPM<IPMDouble> ipm_sdp_solver(constraints.A, constraints.b, constraints.c, &sdp_barrier);
    ipm_sdp_solver.run_solver();
    auto sdp_solution = ipm_sdp_solver.get_solution();
    return ipm_sdp_solver.verify_solution();
}

template<typename IPMDouble>
bool test_lp_solver_random(const int m, const int n) {
    Matrix<IPMDouble> A = Matrix<IPMDouble>::Random(m, n);
    Vector<IPMDouble> x0 = Vector<IPMDouble>::Random(n);
    Vector<IPMDouble> c = Vector<IPMDouble>::Random(n);
    //for loop to create a feasible and bounded instance
    for (int i = 0; i < n; i++) {
        x0(i) += 1;
        c(i) += 1;
    }
    Vector<IPMDouble> b = A * x0;

    LPStandardBarrier<IPMDouble> lp_barrier(n);

    NonSymmetricIPM<IPMDouble> lp_solver(A, b, c, &lp_barrier);
    lp_solver.run_solver();
    Solution<IPMDouble> sol = lp_solver.get_solution();
    return lp_solver.verify_solution(10e-5);
}

template<typename IPMDouble>
void test_sdp_solver(std::string &config_file_str) {

    typedef std::vector<std::pair<IPMDouble, IPMDouble> > HyperRectangle;

    HyperRectangle hyperRectangle;
    hyperRectangle.emplace_back(std::pair<IPMDouble, IPMDouble>(-1, 1));
    const unsigned max_degree = 10;
    const unsigned num_variables = 1;
    EnvelopeProblemSDP <IPMDouble> envelopeProblemSDP(num_variables, max_degree, hyperRectangle);
    Matrix<IPMDouble> poly1 = envelopeProblemSDP.generate_zero_polynomial();
    poly1(0, 0) = 1;
    poly1(1, 0) = 1;
    poly1(-0, 1) = 1;
    poly1(1, 1) = 1;
    poly1(2, 2) = -1;
    poly1(1, 1) = 1;

    Matrix<IPMDouble> poly11 = envelopeProblemSDP.generate_zero_polynomial();
    poly11(0, 0) = 1;
    poly11(1, 0) = -1;
    poly11(-0, 1) = -1;
    poly11(1, 1) = 1;

    envelopeProblemSDP.print_polynomial(poly1);
    envelopeProblemSDP.add_polynomial(poly1);

    envelopeProblemSDP.print_polynomial(poly11);
    envelopeProblemSDP.add_polynomial(poly11);


    Matrix<IPMDouble> p1 = Matrix<IPMDouble>::Zero(poly1.rows(), poly1.cols());
    p1(1, 1) = 10;
    p1(0, 0) = -0;

    Matrix<IPMDouble> p2 = Matrix<IPMDouble>::Zero(poly1.rows(), poly1.cols());
    p2(0, 0) = 10;
    p2(0, 1) = 10;
    for (int i = 2; i <= max_degree; ++i) {
        p2(i, i) = -.0001;
    }
    p2 = p2.eval() + p2.transpose().eval();
    Matrix<IPMDouble> p3 = Matrix<IPMDouble>::Zero(poly1.rows(), poly1.cols());
    p3(0, 0) = 1;
    p3(0, 1) = -1;
    p3 = p3.eval() + p3.transpose().eval();
//
//    envelopeProblemSDP.add_polynomial(p1);
//    envelopeProblemSDP.add_polynomial(p2);
    envelopeProblemSDP.add_polynomial(p3);
//
    Instance<IPMDouble> instance = envelopeProblemSDP.construct_SDP_instance();

//    std::cout << "Objectives Matrix: " << std::endl << envelopeProblemSDP.get_objective_matrix() << std::endl;

    NonSymmetricIPM<IPMDouble> sos_solver(instance, config_file_str);
    sos_solver.run_solver();
    envelopeProblemSDP.print_solution(sos_solver.get_solution());
    envelopeProblemSDP.plot_polynomials_and_solution(sos_solver.get_solution());
}

TEST_CASE ("nonsymmetric_cone_tests") {
    std::cout << "--- Testing LP and SDP solver" << std::endl;
    CHECK(test_lp_solver_random<double>(2, 5));
    CHECK(test_sdp_solver_random_lp_formulation<double>(2, 5));
}
