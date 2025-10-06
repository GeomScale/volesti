// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <list>

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "preprocess/full_dimensional_polytope.hpp"
#include "preprocess/feasible_point.hpp"
#include "convex_bodies/hpolytope.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sampling.hpp"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::SparseMatrix<NT, Eigen::ColMajor> SpMT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

// Helper function to check if a point satisfies Ax <= b
bool satisfies_inequality(const MT& A, const VT& b, const VT& x, NT tol = 1e-6)
{
    VT residual = A * x - b;
    return (residual.array() <= tol).all();
}

// Helper function to check if a point satisfies Aeq*x = beq
bool satisfies_equality(const SpMT& Aeq, const VT& beq, const VT& x, NT tol = 1e-6)
{
    VT residual = Aeq * x - beq;
    return residual.cwiseAbs().maxCoeff() < tol;
}

// Helper function to create sparse matrix from dense
SpMT dense_to_sparse(const MT& dense)
{
    SpMT sparse = dense.sparseView();
    sparse.makeCompressed();
    return sparse;
}

// Helper function to sample from a polytope using AcceleratedBilliardWalk
std::vector<VT> sample_from_polytope(const MT& A, const VT& b, unsigned int num_samples = 100)
{
    unsigned int d = A.cols();
    
    // Create H-polytope
    Hpolytope P(d, A, b);
    
    // Compute inner ball and get its center as starting point
    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();
    Point StartingPoint = inner_ball.first;
    
    // Setup sampling
    RNGType rng(d);
    unsigned int walkL = 1, nburns = 0;
    std::list<Point> randPoints;
    
    // Sample using AcceleratedBilliardWalk
    uniform_sampling<AcceleratedBilliardWalk>(randPoints, P, rng, walkL, num_samples, StartingPoint, nburns);
    
    // Convert to vector of VT
    std::vector<VT> samples;
    samples.reserve(num_samples);
    for (const auto& pt : randPoints)
    {
        samples.push_back(pt.getCoefficients());
    }
    
    return samples;
}

TEST_CASE("full_dimensional_polytope_canonical_simplex")
{
    std::cout << "\n=== Test: Canonical simplex (n-1 dimensional in R^n) ===" << std::endl;
    
    // Canonical simplex in R^3: x1 + x2 + x3 = 1, x1, x2, x3 >= 0
    // This is a 2D simplex embedded in 3D
    int n = 3;
    
    // Equality constraint: x1 + x2 + x3 = 1
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 1;
    
    // Inequality constraints: x1, x2, x3 >= 0
    MT A(n, n);
    A << -1,  0,  0,
          0, -1,  0,
          0,  0, -1;
    
    VT b = VT::Zero(n);
    
    // Compute full dimensional polytope
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    std::cout << "Rank of Aeq: " << Aeq_dense.rows() << std::endl;
    std::cout << "Shift: " << shift.transpose() << std::endl;
    
    // Verify dimension reduction
    CHECK(N.cols() == n - 1);  // Should reduce to 2D
    CHECK(N.rows() == n);
    
    // Verify that N is in the nullspace of Aeq
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    // Verify that the transformation preserves the polytope structure
    CHECK(A_full.rows() == A.rows());
    CHECK(A_full.cols() == N.cols());
    
    // Check that the output polytope is feasible
    VT x_feasible = compute_feasible_point(A_full, b_full);
    CHECK(satisfies_inequality(A_full, b_full, x_feasible, 1e-6));
    std::cout << "Output polytope is feasible" << std::endl;
    
    // Sample from the output (full-dimensional) polytope
    std::cout << "Sampling from output polytope..." << std::endl;
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 50);
    
    // Map sampled points back to original space and verify constraints
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        // Map to original space: y = shift + N * x
        VT y_original = shift + N * x_reduced;
        
        // Check equality constraints
        bool eq_satisfied = satisfies_equality(Aeq, beq, y_original, 1e-5);
        
        // Check inequality constraints  
        bool ineq_satisfied = satisfies_inequality(A, b, y_original, 1e-5);
        
        if (eq_satisfied && ineq_satisfied)
            num_valid++;
    }
    
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_simple_feasible")
{
    std::cout << "\n=== Test: Simple feasible case (2D in 3D) ===" << std::endl;
    
    // 2D plane in 3D: x + y + z = 3, with bounds 0 <= x,y,z <= 2
    int n = 3;
    
    // Equality constraint: x + y + z = 3
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 3;
    
    // Inequality constraints: 0 <= x,y,z <= 2
    MT A(6, n);
    A << -1,  0,  0,
          0, -1,  0,
          0,  0, -1,
          1,  0,  0,
          0,  1,  0,
          0,  0,  1;
    
    VT b(6);
    b << 0, 0, 0, 2, 2, 2;
    
    // Compute full dimensional polytope
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    std::cout << "Shift: " << shift.transpose() << std::endl;
    
    // Verify dimension reduction
    CHECK(N.cols() == n - 1);
    
    // Verify N is in nullspace
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    // Check output polytope feasibility
    VT x_feasible = compute_feasible_point(A_full, b_full);
    CHECK(satisfies_inequality(A_full, b_full, x_feasible, 1e-6));
    
    // Sample and verify transformation
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 50);
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        VT y_original = shift + N * x_reduced;
        if (satisfies_equality(Aeq, beq, y_original, 1e-5) && 
            satisfies_inequality(A, b, y_original, 1e-5))
            num_valid++;
    }
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_line_segment")
{
    std::cout << "\n=== Test: Line segment (1D in 2D) ===" << std::endl;
    
    int n = 2;
    
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 1;
    
    MT A(4, n);
    A << -1,  0,
          1,  0,
          0, -1,
          0,  1;
    
    VT b(4);
    b << 0, 1, 0, 1;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    
    CHECK(N.cols() == 1);
    
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    VT x_feasible = compute_feasible_point(A_full, b_full);
    CHECK(satisfies_inequality(A_full, b_full, x_feasible, 1e-6));
    
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 50);
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        VT y_original = shift + N * x_reduced;
        if (satisfies_equality(Aeq, beq, y_original, 1e-5) && 
            satisfies_inequality(A, b, y_original, 1e-5))
            num_valid++;
    }
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_high_dimensional")
{
    std::cout << "\n=== Test: High dimensional (8D in 10D) ===" << std::endl;
    
    int n = 10;
    int m = 2;
    
    MT Aeq_dense(m, n);
    Aeq_dense << 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 1, 1, 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(m);
    beq << 5, 5;
    
    MT A = -MT::Identity(n, n);
    VT b = VT::Zero(n);
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    
    CHECK(N.cols() == n - m);
    
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    VT x_feasible = compute_feasible_point(A_full, b_full);
    CHECK(satisfies_inequality(A_full, b_full, x_feasible, 1e-6));
    
    // For high dimensional, test with fewer samples
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 30);
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        VT y_original = shift + N * x_reduced;
        if (satisfies_equality(Aeq, beq, y_original, 1e-5) && 
            satisfies_inequality(A, b, y_original, 1e-5))
            num_valid++;
    }
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_multiple_constraints")
{
    std::cout << "\n=== Test: Multiple equality constraints (3D to 1D) ===" << std::endl;
    
    int n = 3;
    int m = 2;
    
    MT Aeq_dense(m, n);
    Aeq_dense << 1, 1, 0,
                 0, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(m);
    beq << 1, 1;
    
    MT A(6, n);
    A << -1,  0,  0,
          0, -1,  0,
          0,  0, -1,
          1,  0,  0,
          0,  1,  0,
          0,  0,  1;
    
    VT b(6);
    b << 0, 0, 0, 1, 1, 1;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    
    CHECK(N.cols() == n - m);
    CHECK(N.cols() == 1);
    
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    VT x_feasible = compute_feasible_point(A_full, b_full);
    CHECK(satisfies_inequality(A_full, b_full, x_feasible, 1e-6));
    
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 50);
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        VT y_original = shift + N * x_reduced;
        if (satisfies_equality(Aeq, beq, y_original, 1e-5) && 
            satisfies_inequality(A, b, y_original, 1e-5))
            num_valid++;
    }
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_hypercube_intersection")
{
    std::cout << "\n=== Test: Hypercube with plane (3D to 2D) ===" << std::endl;
    
    int n = 3;
    
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 1.5;
    
    MT A(6, n);
    A << -1,  0,  0,
          0, -1,  0,
          0,  0, -1,
          1,  0,  0,
          0,  1,  0,
          0,  0,  1;
    
    VT b(6);
    b << 0, 0, 0, 1, 1, 1;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    
    CHECK(N.cols() == 2);
    
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    VT x_feasible = compute_feasible_point(A_full, b_full);
    CHECK(satisfies_inequality(A_full, b_full, x_feasible, 1e-6));
    
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 50);
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        VT y_original = shift + N * x_reduced;
        if (satisfies_equality(Aeq, beq, y_original, 1e-5) && 
            satisfies_inequality(A, b, y_original, 1e-5))
            num_valid++;
    }
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_sparse_constraint")
{
    std::cout << "\n=== Test: Sparse equality constraint (5D to 4D) ===" << std::endl;
    
    int n = 5;
    
    MT Aeq_dense = MT::Zero(1, n);
    Aeq_dense(0, 0) = 1;
    Aeq_dense(0, 4) = 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 2;
    
    MT A(2*n, n);
    A.topRows(n) = -MT::Identity(n, n);
    A.bottomRows(n) = MT::Identity(n, n);
    
    VT b(2*n);
    b.head(n).setZero();
    b.tail(n).setConstant(2);
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    
    CHECK(N.cols() == 4);
    
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    VT x_feasible = compute_feasible_point(A_full, b_full);
    CHECK(satisfies_inequality(A_full, b_full, x_feasible, 1e-6));
    
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 50);
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        VT y_original = shift + N * x_reduced;
        if (satisfies_equality(Aeq, beq, y_original, 1e-5) && 
            satisfies_inequality(A, b, y_original, 1e-5))
            num_valid++;
    }
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_orthogonality")
{
    std::cout << "\n=== Test: Nullspace orthogonality properties ===" << std::endl;
    
    int n = 4;
    
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 4;
    
    MT A = -MT::Identity(n, n);
    VT b = VT::Zero(n);
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    std::cout << "Original dimension: " << n << std::endl;
    std::cout << "Reduced dimension: " << N.cols() << std::endl;
    
    MT Aeq_N = Aeq_dense * N;
    CHECK(Aeq_N.cwiseAbs().maxCoeff() < 1e-6);
    
    MT NtN = N.transpose() * N;
    MT I = MT::Identity(N.cols(), N.cols());
    CHECK((NtN - I).cwiseAbs().maxCoeff() < 1e-6);
    
    std::cout << "N^T * N =\n" << NtN << std::endl;
    std::cout << "Columns are orthonormal: " 
              << ((NtN - I).cwiseAbs().maxCoeff() < 1e-6 ? "YES" : "NO") << std::endl;
    
    // Also verify with sampling
    std::vector<VT> samples = sample_from_polytope(A_full, b_full, 30);
    int num_valid = 0;
    for (const auto& x_reduced : samples)
    {
        VT y_original = shift + N * x_reduced;
        if (satisfies_equality(Aeq, beq, y_original, 1e-5) && 
            satisfies_inequality(A, b, y_original, 1e-5))
            num_valid++;
    }
    std::cout << "Valid samples: " << num_valid << " / " << samples.size() << std::endl;
    CHECK(num_valid == samples.size());
}

TEST_CASE("full_dimensional_polytope_infeasible")
{
    std::cout << "\n=== Test: Infeasible system (plane outside cube) ===" << std::endl;
    
    int n = 3;
    
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 5;  // Impossible since max(x+y+z) = 3 in unit cube
    
    MT A(6, n);
    A << -1,  0,  0,
          0, -1,  0,
          0,  0, -1,
          1,  0,  0,
          0,  1,  0,
          0,  0,  1;
    
    VT b(6);
    b << 0, 0, 0, 1, 1, 1;
    
    // This should still produce output, but the output polytope should be infeasible
    try {
        auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
        
        std::cout << "Shift found: " << shift.transpose() << std::endl;
        
        // Check if equality is satisfied by shift (it should be)
        CHECK(satisfies_equality(Aeq, beq, shift, 1e-6));
        std::cout << "Equality constraint satisfied by shift: YES" << std::endl;
        
        // Shift should NOT satisfy the inequalities (problem is infeasible)
        bool shift_feasible = satisfies_inequality(A, b, shift, 1e-6);
        std::cout << "Shift satisfies inequalities: " << (shift_feasible ? "YES" : "NO") << std::endl;
        CHECK(!shift_feasible);  // Should be infeasible
        
        // The output polytope A_full * x <= b_full should also be infeasible
        std::cout << "Attempting to find feasible point in output polytope..." << std::endl;
        bool output_polytope_feasible = true;
        try {
            VT x_feasible = compute_feasible_point(A_full, b_full);
            // If we reach here, check if it actually satisfies the constraints
            if (satisfies_inequality(A_full, b_full, x_feasible, 1e-6)) {
                std::cout << "WARNING: Output polytope appears feasible!" << std::endl;
                output_polytope_feasible = true;
            } else {
                output_polytope_feasible = false;
            }
        } catch (const std::exception& e) {
            std::cout << "compute_feasible_point failed (expected): " << e.what() << std::endl;
            output_polytope_feasible = false;
        }
        
        // For infeasible input, output polytope should also be infeasible
        std::cout << "Output polytope feasible: " << (output_polytope_feasible ? "YES" : "NO") << std::endl;
        CHECK(!output_polytope_feasible);
        
    } catch (const std::exception& e) {
        std::cout << "Exception during transformation: " << e.what() << std::endl;
    }
}
