// VolEsti (volume computation and sampling library)
// Licensed under GNU LGPL.3, see LICENCE file

// Tests for OrderPolytope trigonometric and quadratic boundary oracles.
// These tests validate the new HMC boundary oracles added in Phase 1
// without requiring any volume computation (avoiding MSVC khach.h issue).

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "doctest.h"
#include <iostream>
#include <cmath>

#include <boost/random.hpp>

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "cartesian_geom/point.h"
#include "convex_bodies/orderpolytope.h"
#include "misc/poset.h"


// ------------------------------------------------------------------
// Test: trigonometric_positive_intersect for OrderPolytope
// Validates that the Gaussian HMC boundary oracle finds a valid
// intersection time with the polytope boundary.
// ------------------------------------------------------------------
template <typename NT>
void call_test_trigonometric_intersect() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef typename OrderPolytope<Point>::VT VT;
    typedef typename Poset::RV RV;

    // Create Poset: 4 elements, a0 <= a1, a0 <= a2, a1 <= a3
    RV poset_data{{0, 1}, {0, 2}, {1, 3}};
    Poset poset(4, poset_data);

    OrderPolytope<Point> OP(poset);
    unsigned int d = OP.dimension();

    // Start from interior point
    VT ip = OP.inner_point();
    Point r(ip);

    // Direction
    Point v(d);
    for (unsigned int i = 0; i < d; i++)
        v.set_coord(i, NT(i + 1) * 0.3);

    NT omega = 1.0;
    int facet_prev = -1;

    auto result = OP.trigonometric_positive_intersect(r, v, omega, facet_prev);
    NT t_hit = result.first;
    int facet_hit = result.second;

    // t should be positive and finite
    CHECK(t_hit > NT(0));
    CHECK(t_hit < NT(1e15));

    // facet should be valid
    CHECK(facet_hit >= 0);
    CHECK(facet_hit < (int)OP.num_of_hyperplanes());

    // The point at t_hit should be on the boundary
    NT sinVal = std::sin(omega * t_hit);
    NT cosVal = std::cos(omega * t_hit);
    Point p_boundary(d);
    for (unsigned int i = 0; i < d; i++) {
        p_boundary.set_coord(i, cosVal * r[i] + (sinVal / omega) * v[i]);
    }

    // Check the boundary point is inside or on the boundary (with tolerance)
    int in_result = OP.is_in(p_boundary, NT(1e-6));
    CHECK((in_result == -1 || in_result == 0));

    std::cout << "Trigonometric intersect: t=" << t_hit
              << ", facet=" << facet_hit << std::endl;
}


// ------------------------------------------------------------------
// Test: quadratic_positive_intersect for OrderPolytope
// Validates that the Exponential HMC boundary oracle works correctly.
// ------------------------------------------------------------------
template <typename NT>
void call_test_quadratic_intersect() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef typename OrderPolytope<Point>::VT VT;
    typedef typename Poset::RV RV;

    // Chain: a0 <= a1 <= a2
    RV poset_data{{0, 1}, {1, 2}};
    Poset poset(3, poset_data);

    OrderPolytope<Point> OP(poset);
    unsigned int d = OP.dimension();
    unsigned int m = OP.num_of_hyperplanes();

    VT ip = OP.inner_point();
    Point r(ip);
    Point v = Point::all_ones(d);

    // Bias vector and temperature for exponential distribution
    Point c = Point::all_ones(d);
    NT T = 1.0;

    VT Ac = OP.vec_mult(c.getCoefficients());
    VT Ar = VT::Zero(m);
    VT Av = VT::Zero(m);
    int facet_prev = -1;

    auto result = OP.quadratic_positive_intersect(r, v, Ac, T, Ar, Av, facet_prev);
    NT t_hit = result.first;
    int facet_hit = result.second;

    CHECK(t_hit > NT(0));
    CHECK(t_hit < NT(1e15));
    CHECK(facet_hit >= 0);
    CHECK(facet_hit < (int)m);

    std::cout << "Quadratic intersect: t=" << t_hit
              << ", facet=" << facet_hit << std::endl;
}


// ------------------------------------------------------------------
// Test: vec_mult consistency (A*x via sparse vs expected)
// Verifies that the OrderPolytope's vec_mult gives correct results
// for the new oracle methods.
// ------------------------------------------------------------------
template <typename NT>
void call_test_vec_mult_for_oracles() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef typename OrderPolytope<Point>::VT VT;
    typedef typename OrderPolytope<Point>::MT MT;
    typedef typename Poset::RV RV;

    // b < a, c   (a=0, b=1, c=2, so 1<=0, 1<=2)
    RV poset_data{{1, 0}, {1, 2}};
    Poset poset(3, poset_data);

    OrderPolytope<Point> OP(poset);
    unsigned int d = OP.dimension();
    unsigned int m = OP.num_of_hyperplanes();

    // Test point
    VT x(d);
    x << 0.7, 0.3, 0.8;

    // Compute A*x using vec_mult (optimized sparse)
    VT Ax_sparse = OP.vec_mult(x);

    // Compute A*x using dense multiplication
    MT A = OP.get_dense_mat();
    VT Ax_dense = A * x;

    // They should be equal
    CHECK(Ax_sparse.size() == Ax_dense.size());
    for (unsigned int i = 0; i < (unsigned int)Ax_sparse.size(); i++) {
        CHECK(std::abs(Ax_sparse(i) - Ax_dense(i)) < NT(1e-10));
    }

    std::cout << "vec_mult consistency: PASSED" << std::endl;
}


// ------------------------------------------------------------------
// Test: Oracle gives consistent results for different poset types
// ------------------------------------------------------------------
template <typename NT>
void call_test_oracle_antichain() {
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef typename OrderPolytope<Point>::VT VT;
    typedef typename Poset::RV RV;

    // Antichain: 3 elements, no relations
    // OrderPolytope = [0,1]^3 (unit cube)
    RV poset_data{};
    Poset poset(3, poset_data);

    OrderPolytope<Point> OP(poset);
    unsigned int d = OP.dimension();

    // Point at center of cube
    Point r(d);
    r.set_coord(0, 0.5);
    r.set_coord(1, 0.5);
    r.set_coord(2, 0.5);

    // Direction towards (1,1,1)
    Point v(d);
    v.set_coord(0, 1.0);
    v.set_coord(1, 1.0);
    v.set_coord(2, 1.0);

    NT omega = 1.0;
    int facet_prev = -1;

    auto result = OP.trigonometric_positive_intersect(r, v, omega, facet_prev);
    NT t_hit = result.first;
    int facet_hit = result.second;

    CHECK(t_hit > NT(0));
    CHECK(facet_hit >= 0);

    std::cout << "Antichain oracle: t=" << t_hit
              << ", facet=" << facet_hit << std::endl;
}


// ------------------------------------------------------------------
// Register test cases
// ------------------------------------------------------------------
TEST_CASE("trigonometric_intersect") {
    call_test_trigonometric_intersect<double>();
}

TEST_CASE("quadratic_intersect") {
    call_test_quadratic_intersect<double>();
}

TEST_CASE("vec_mult_for_oracles") {
    call_test_vec_mult_for_oracles<double>();
}

TEST_CASE("oracle_antichain") {
    call_test_oracle_antichain<double>();
}
