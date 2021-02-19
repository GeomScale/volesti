// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "sos/utils.h"
#include "sos/NonSymmetricIPM.h"
#include "sos/barriers/LPStandardBarrier.h"

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

TEST_CASE ("nonsymmetric_lp_tests") {
    std::cout << "--- Test Linear Programs with Nonsymmetric Conic Solver" << std::endl;
    CHECK(test_lp_solver_random<double>(2, 5));
}
