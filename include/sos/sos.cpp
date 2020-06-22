#include "NonSymmetricIPM.h"
#include "MonomialsClass.h"
#include "EnvelopeProblemSDP.h"
#include "EnvelopeProblemSOS.h"
#include "spdlog/spdlog.h"
#include "spdlog/cfg/env.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include <fstream>

//TODO: Use ARPACK for faster sparse linear algebra.

Constraints convert_LP_to_SDP(Matrix & A, Vector & b, Vector & c) {
    const int m = A.rows();
    const int n = A.cols();

    Constraints SDP_constraints;
    SDP_constraints.A = Matrix::Zero(m + n * n - n, n * n);

    std::cout << A << std::endl;
    for (int i = 0; i < m; ++i) {
        Matrix diag_row = A.row(i).asDiagonal();
        Eigen::Map<Matrix> diag_row_stacked(diag_row.data(), 1, n * n);
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

    SDP_constraints.b = Vector::Zero(SDP_constraints.A.rows());
    SDP_constraints.b.block(0, 0, m, 1) = b;

    SDP_constraints.c = Vector::Zero(SDP_constraints.A.cols());
    for (int j = 0; j < n; ++j) {
        SDP_constraints.c(j + j * n) = c(j);
    }

    return SDP_constraints;
}

bool test_sdp_solver_random_lp_formulation(const int m, const int n){
    Matrix A = Matrix::Random(m, n);
    Vector x0 = Vector::Random(n);
    Vector c = Vector::Random(n);
    //for loop to create a feasible and bounded instance
    for (int i = 0; i < n; i++) {
        x0(i) += 1;
        c(i) += 1;
    }
    Vector b = A * x0;

    Constraints constraints = convert_LP_to_SDP(A, b, c);

    SDPStandardBarrier sdp_barrier(n);

    NonSymmetricIPM ipm_sdp_solver(constraints.A, constraints.b, constraints.c, &sdp_barrier);
    std::clock_t sdp_timer = std::clock();
    ipm_sdp_solver.run_solver();
    auto sdp_solution = ipm_sdp_solver.get_solution();
    return ipm_sdp_solver.verify_solution();
}

bool test_lp_solver_random(const int m, const int n){
    Matrix A = Matrix::Random(m, n);
    Vector x0 = Vector::Random(n);
    Vector c = Vector::Random(n);
    //for loop to create a feasible and bounded instance
    for (int i = 0; i < n; i++) {
        x0(i) += 1;
        c(i) += 1;
    }
    Vector b = A * x0;

    LPStandardBarrier lp_barrier(n);

    NonSymmetricIPM lp_solver(A, b, c, &lp_barrier);
    lp_solver.run_solver();
    Solution sol = lp_solver.get_solution();
    return lp_solver.verify_solution(10e-5);
}

int main(int const argc, char **argv) {
    std::ifstream file;
    if(argc < 2){
        std::cout << "No file provided. Use default file instead." << std::endl;
        file.open("../config/default.txt");
        if(not file.is_open()){
            file.open("config/default.txt");
        } if(not file.is_open()){
            std::cout << "Could not locate file." << std::endl;
            return 1;
        }

    } else{
        file.open(argv[1]);
    }

    srand(time(nullptr));

    assert(test_lp_solver_random(2,5));
    assert(test_sdp_solver_random_lp_formulation(2,5));

    spdlog::logger logger("main_logger");

    logger.set_level(spdlog::level::info);
    logger.info("Run Solvers!");

    // create color multi threaded logger
    auto console = spdlog::stdout_color_mt("console");
    console->info("Logger level is {}", console->level());

    std::string line;
    std::getline(file, line);
    std::istringstream  iss(line);
    int max_degree;
    iss >> max_degree;

    HyperRectangle hyperRectangle;
    //Note: Keep interval bounds
    hyperRectangle.push_back(std::pair<Double, Double>(-1, 1));
//    EnvelopeProblemSDP envelopeProblem(1, max_degree, hyperRectangle);
    EnvelopeProblemSOS envelopeProblemSos(1, max_degree, hyperRectangle);

    while(std::getline(file, line)){
        std::cout << "Read line: " << line << std::endl;
        std::istringstream poly_stream(line);
        PolynomialSOS sos_poly = envelopeProblemSos.generate_zero_polynomial();
        Double val;
        unsigned idx = 0;
        while(poly_stream >> val){
            if(idx >= sos_poly.size()){
                return 1;
            }
            sos_poly[idx++] = val;
        }
        std::cout << sos_poly.transpose() << std::endl;
        envelopeProblemSos.add_polynomial(sos_poly);
    }

    Instance instance_interp = envelopeProblemSos.construct_SOS_instance();

    NonSymmetricIPM sos_solver_interp(instance_interp);

    sos_solver_interp.run_solver();

    envelopeProblemSos.print_solution(sos_solver_interp.get_solution());
    envelopeProblemSos.plot_polynomials_and_solution(sos_solver_interp.get_solution());
    return 0;

//    PolynomialSDP poly1 = envelopeProblem.generate_zero_polynomial();
//    poly1(0,0) = 1;
//    poly1(1, 0) = 1;
//    poly1(-0,1) = 1;
//    poly1(1, 1) = 1;
//    poly1(2,2) = -1;
//    poly1(1,1) = 1;
//
//    Polynomial poly11 = envelopeProblem.generate_zero_polynomial();
//    poly11(0,0) = 1;
//    poly11(1, 0) = -1;
//    poly11(-0,1) = -1;
//    poly11(1, 1) = 1;
//
//    envelopeProblem.print_polynomial(poly1);
//    envelopeProblem.add_polynomial(poly1);
//
//    envelopeProblem.print_polynomial(poly11);
//    envelopeProblem.add_polynomial(poly11);


//    PolynomialSDP p1 = Matrix::Zero(poly1.rows(), poly1.cols());
//    p1(1,1) = 10;
//    p1(0, 0) = -0;
//
//    PolynomialSDP p2 = Matrix::Zero(poly1.rows(), poly1.cols());
//    p2(0,0) = 10;
//    p2(0,1) = 10;
//    for (int i = 2; i <= max_degree; ++i) {
//        p2(i,i) = -.0001;
//    }
//    p2 = p2.eval() + p2.transpose().eval();
//    PolynomialSDP p3 = Matrix::Zero(poly1.rows(), poly1.cols());
//    p3(0,0) = 1;
//    p3(0,1) = -1;
//    p3 = p3.eval() + p3.transpose().eval();
////
//    envelopeProblem.add_polynomial(p1);
//    envelopeProblem.add_polynomial(p2);
//    envelopeProblem.add_polynomial(p3);
////
//    Instance instance = envelopeProblem.construct_SDP_instance();
//
//    std::cout << "Objectives Matrix: " << std::endl << envelopeProblem.get_objective_matrix() << std::endl;
////
//    NonSymmetricIPM sos_solver(instance);
//    sos_solver.run_solver();
//    envelopeProblem.print_solution(sos_solver.get_solution());
//    envelopeProblem.plot_polynomials_and_solution(sos_solver.get_solution());
//    return 0;

}
