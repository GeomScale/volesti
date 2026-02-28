// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include "preprocess/full_dimensional_polytope.hpp"
#include "generators/known_polytope_generators.h"
#include "convex_bodies/hpolytope.h"
#include "cartesian_geom/cartesian_kernel.h"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::SparseMatrix<NT, Eigen::ColMajor> SpMT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;

// Helper function to create sparse matrix from dense
SpMT dense_to_sparse(const MT& dense)
{
    SpMT sparse = dense.sparseView();
    sparse.makeCompressed();
    return sparse;
}

// Helper function to print section header
void print_header(const std::string& title)
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "  " << title << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}

// Helper function to print results
void print_results(int original_dim, int reduced_dim, const VT& shift, const MT& N, 
                   const MT& A_full, const VT& b_full)
{
    std::cout << "\nRESULTS:" << std::endl;
    std::cout << "--------" << std::endl;
    std::cout << "Original dimension: " << original_dim << std::endl;
    std::cout << "Reduced dimension:  " << reduced_dim << std::endl;
    std::cout << "Dimension reduction: " << original_dim << "D -> " << reduced_dim << "D" << std::endl;
    
    std::cout << "\nShift vector (point satisfying Aeq*x = beq):" << std::endl;
    std::cout << shift.transpose() << std::endl;
    
    std::cout << "\nNullspace basis N (" << N.rows() << " x " << N.cols() << "):" << std::endl;
    std::cout << N << std::endl;
    
    std::cout << "\nOutput polytope: A_full * y <= b_full" << std::endl;
    std::cout << "  A_full dimensions: " << A_full.rows() << " x " << A_full.cols() << std::endl;
    std::cout << "  Number of constraints: " << A_full.rows() << std::endl;
    
    // Verify orthonormality
    MT NtN = N.transpose() * N;
    MT I = MT::Identity(N.cols(), N.cols());
    NT orthogonality_error = (NtN - I).cwiseAbs().maxCoeff();
    std::cout << "\nNullspace orthonormality check:" << std::endl;
    std::cout << "  ||N^T * N - I||_inf = " << orthogonality_error << std::endl;
    std::cout << "  Columns are orthonormal: " << (orthogonality_error < 1e-6 ? "YES" : "NO") << std::endl;
}

/**
 * Example 1: Canonical Simplex
 * The standard (n-1)-simplex embedded in n-dimensional space
 * Equality: x1 + x2 + ... + xn = 1
 * Inequalities: xi >= 0 for all i
 */
void example_canonical_simplex()
{
    print_header("Example 1: Canonical Simplex (2D triangle in 3D space)");
    
    std::cout << "\nDESCRIPTION:" << std::endl;
    std::cout << "The canonical 2-simplex (triangle) in 3D satisfies:" << std::endl;
    std::cout << "  Equality constraint:   x + y + z = 1" << std::endl;
    std::cout << "  Inequality constraints: x, y, z >= 0" << std::endl;
    std::cout << "This is a 2D surface (triangle) embedded in 3D space." << std::endl;
    
    int n = 3;
    
    // Equality constraint: x + y + z = 1
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 1;
    
    // Inequality constraints: x, y, z >= 0
    MT A(n, n);
    A << -1,  0,  0,
          0, -1,  0,
          0,  0, -1;
    
    VT b = VT::Zero(n);
    
    std::cout << "\nINPUT:" << std::endl;
    std::cout << "  Aeq (equality constraints):" << std::endl;
    std::cout << "    " << Aeq_dense << std::endl;
    std::cout << "  beq:" << std::endl;
    std::cout << "    " << beq.transpose() << std::endl;
    
    // Compute full dimensional polytope
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    print_results(n, N.cols(), shift, N, A_full, b_full);
}

/**
 * Example 2: Hyperplane intersecting a cube
 * A plane cutting through a unit cube
 */
void example_plane_cube_intersection()
{
    print_header("Example 2: Plane intersecting unit cube");
    
    std::cout << "\nDESCRIPTION:" << std::endl;
    std::cout << "A plane x + y + z = 1.5 intersecting the unit cube [0,1]^3." << std::endl;
    std::cout << "The intersection is a 2D hexagon embedded in 3D." << std::endl;
    
    int n = 3;
    
    // Equality constraint: x + y + z = 1.5
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 1.5;
    
    // Inequality constraints: 0 <= x, y, z <= 1
    MT A(6, n);
    A << -1,  0,  0,
          0, -1,  0,
          0,  0, -1,
          1,  0,  0,
          0,  1,  0,
          0,  0,  1;
    
    VT b(6);
    b << 0, 0, 0, 1, 1, 1;
    
    std::cout << "\nINPUT:" << std::endl;
    std::cout << "  Equality: x + y + z = 1.5" << std::endl;
    std::cout << "  Inequalities: 0 <= x, y, z <= 1 (unit cube)" << std::endl;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    print_results(n, N.cols(), shift, N, A_full, b_full);
}

/**
 * Example 3: Multiple equality constraints
 * Reducing 4D to 2D with two independent constraints
 */
void example_multiple_constraints()
{
    print_header("Example 3: Multiple equality constraints (4D -> 2D)");
    
    std::cout << "\nDESCRIPTION:" << std::endl;
    std::cout << "Two independent equality constraints in 4D:" << std::endl;
    std::cout << "  x1 + x2 = 2" << std::endl;
    std::cout << "  x3 + x4 = 3" << std::endl;
    std::cout << "With bounds: 0 <= xi <= 2 for all i" << std::endl;
    std::cout << "This reduces the problem from 4D to 2D." << std::endl;
    
    int n = 4;
    int m = 2;
    
    // Two equality constraints
    MT Aeq_dense(m, n);
    Aeq_dense << 1, 1, 0, 0,
                 0, 0, 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(m);
    beq << 2, 3;
    
    // Inequality constraints: 0 <= xi <= 2
    MT A(2*n, n);
    A.topRows(n) = -MT::Identity(n, n);
    A.bottomRows(n) = MT::Identity(n, n);
    
    VT b(2*n);
    b.head(n).setZero();
    b.tail(n).setConstant(2);
    
    std::cout << "\nINPUT:" << std::endl;
    std::cout << "  Aeq (two equality constraints):" << std::endl;
    std::cout << Aeq_dense << std::endl;
    std::cout << "  beq: " << beq.transpose() << std::endl;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    print_results(n, N.cols(), shift, N, A_full, b_full);
}

/**
 * Example 4: Sparse equality constraint
 * Demonstrating efficient handling of sparse constraints
 */
void example_sparse_constraint()
{
    print_header("Example 4: Sparse equality constraint (5D -> 4D)");
    
    std::cout << "\nDESCRIPTION:" << std::endl;
    std::cout << "A sparse constraint involving only first and last variables:" << std::endl;
    std::cout << "  x1 + x5 = 2" << std::endl;
    std::cout << "With bounds: 0 <= xi <= 2 for all i" << std::endl;
    std::cout << "This is common in network flow or coupling constraints." << std::endl;
    
    int n = 5;
    
    // Sparse equality: only x1 and x5
    MT Aeq_dense = MT::Zero(1, n);
    Aeq_dense(0, 0) = 1;
    Aeq_dense(0, 4) = 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 2;
    
    // Box constraints
    MT A(2*n, n);
    A.topRows(n) = -MT::Identity(n, n);
    A.bottomRows(n) = MT::Identity(n, n);
    
    VT b(2*n);
    b.head(n).setZero();
    b.tail(n).setConstant(2);
    
    std::cout << "\nINPUT:" << std::endl;
    std::cout << "  Sparse equality: x1 + x5 = 2" << std::endl;
    std::cout << "  (variables x2, x3, x4 are unconstrained by equality)" << std::endl;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    print_results(n, N.cols(), shift, N, A_full, b_full);
}

/**
 * Example 5: Line segment in 2D
 * Extreme case: reducing to 1D
 */
void example_line_segment()
{
    print_header("Example 5: Line segment (1D line in 2D plane)");
    
    std::cout << "\nDESCRIPTION:" << std::endl;
    std::cout << "The line x + y = 1 within the square [0,1]^2." << std::endl;
    std::cout << "This is the most extreme reduction: 2D -> 1D." << std::endl;
    
    int n = 2;
    
    // Equality: x + y = 1
    MT Aeq_dense(1, n);
    Aeq_dense << 1, 1;
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 1;
    
    // Box constraints
    MT A(4, n);
    A << -1,  0,
          1,  0,
          0, -1,
          0,  1;
    
    VT b(4);
    b << 0, 1, 0, 1;
    
    std::cout << "\nINPUT:" << std::endl;
    std::cout << "  Equality: x + y = 1" << std::endl;
    std::cout << "  Inequalities: 0 <= x, y <= 1" << std::endl;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    print_results(n, N.cols(), shift, N, A_full, b_full);
}

/**
 * Example 6: Birkhoff polytope
 * Real-world example: doubly stochastic matrices
 */
void example_birkhoff_polytope()
{
    print_header("Example 6: Birkhoff Polytope (doubly stochastic matrices)");
    
    std::cout << "\nDESCRIPTION:" << std::endl;
    std::cout << "The Birkhoff polytope B_3 consists of 3x3 doubly stochastic matrices:" << std::endl;
    std::cout << "  - Each row sums to 1" << std::endl;
    std::cout << "  - Each column sums to 1" << std::endl;
    std::cout << "  - All entries are non-negative" << std::endl;
    std::cout << "This is a 9D space with 6 equality constraints -> 3D polytope" << std::endl;
    
    // Generate Birkhoff polytope of size 3
    Hpolytope birkhoff = generate_birkhoff<Hpolytope>(3);
    
    int n = birkhoff.dimension();  // Should be 9 for 3x3 matrices
    MT A = birkhoff.get_mat();
    VT b = birkhoff.get_vec();
    
    std::cout << "\nOriginal polytope dimension: " << n << std::endl;
    std::cout << "Number of inequality constraints: " << A.rows() << std::endl;
    
    // For demonstration, let's use a simple equality: sum of all entries = 3
    MT Aeq_dense = MT::Ones(1, n);
    SpMT Aeq = dense_to_sparse(Aeq_dense);
    
    VT beq(1);
    beq << 3;  // Sum of all entries in a 3x3 doubly stochastic matrix
    
    std::cout << "\nUsing simple equality constraint: sum of all entries = 3" << std::endl;
    
    auto [A_full, b_full, shift, N] = compute_full_dimensional_polytope(Aeq, beq, A, b);
    
    print_results(n, N.cols(), shift, N, A_full, b_full);
}

int main()
{
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "  FULL DIMENSIONAL POLYTOPE EXAMPLES" << std::endl;
    std::cout << "  Transforming lower-dimensional polytopes to full dimension" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nOVERVIEW:" << std::endl;
    std::cout << "These examples demonstrate the compute_full_dimensional_polytope function," << std::endl;
    std::cout << "which transforms a polytope defined by:" << std::endl;
    std::cout << "  - Equality constraints: Aeq * x = beq" << std::endl;
    std::cout << "  - Inequality constraints: A * x <= b" << std::endl;
    std::cout << "into a full-dimensional polytope: A_full * y <= b_full" << std::endl;
    std::cout << "\nThe transformation provides:" << std::endl;
    std::cout << "  - shift: A point satisfying the equality constraints" << std::endl;
    std::cout << "  - N: Nullspace basis (orthonormal columns)" << std::endl;
    std::cout << "  - Mapping: x = shift + N * y" << std::endl;
    
    // Run all examples
    example_canonical_simplex();
    std::this_thread::sleep_for(std::chrono::seconds(1));

    example_plane_cube_intersection();
    std::this_thread::sleep_for(std::chrono::seconds(1));

    example_multiple_constraints();
    std::this_thread::sleep_for(std::chrono::seconds(1));

    example_sparse_constraint();
    std::this_thread::sleep_for(std::chrono::seconds(1));

    example_line_segment();
    std::this_thread::sleep_for(std::chrono::seconds(1));

    example_birkhoff_polytope();
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "  ALL EXAMPLES COMPLETED" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    return 0;
}
