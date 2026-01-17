#include "doctest.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>        // ✅ REQUIRED
#include <Eigen/SparseCholesky> // ✅ REQUIRED

#include <cmath>

#include "preprocess/max_inscribed_ball.hpp"

TEST_CASE("max_inscribed_ball: unit box [-1,1]^n")
{
    using NT = double;
    using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
    using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

    const int n = 5;
    const NT tol = 1e-8;

    MT A(2 * n, n);
    VT b(2 * n);

    A.topRows(n) = MT::Identity(n, n);
    A.bottomRows(n) = -MT::Identity(n, n);
    b.setOnes();

    auto [x, r, converged] = max_inscribed_ball(A, b, 10000, tol);

    CHECK(converged);
    CHECK(std::abs(r - 1.0) < 1e-6);
    CHECK(x.norm() < 1e-6);
    CHECK((A * x - b).maxCoeff() <= 1e-6);
}

TEST_CASE("max_inscribed_ball invalidates radius on non-convergence") {
    using NT = double;
    using MT = Eigen::MatrixXd;
    using VT = Eigen::VectorXd;

    // Define a very thin box in 2D
    MT A(4, 2);
    VT b(4);

    A <<  1,  0,
         -1,  0,
          0,  1,
          0, -1;

    // Extremely thin in x-direction → likely non-convergence
    b << 1e-12, 1e-12, 1.0, 1.0;

    // Force non-convergence
    unsigned int maxiter = 1;
    NT tol = 1e-12;

    auto [x, t, converge] = max_inscribed_ball(A, b, maxiter, tol);

    REQUIRE_FALSE(converge);
    REQUIRE(t < NT(0));   // radius must be invalidated
}