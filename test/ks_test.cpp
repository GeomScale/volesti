// test/ks_test.cpp
#include "doctest.h"
#include <random>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "diagnostics/KS_test.hpp"

typedef Cartesian<double>                  Kernel;
typedef typename Kernel::Point             Point;
typedef HPolytope<Point>                   HPOLYTOPE;
typedef typename HPOLYTOPE::MT             MT;
typedef typename HPOLYTOPE::VT             VT;

static HPOLYTOPE make_simplex(int d) {
    MT A = MT::Zero(d + 1, d);
    VT b = VT::Zero(d + 1);
    for (int i = 0; i < d; ++i) A(i, i) = -1.0;
    A.row(d).setOnes();                          
    b(d) = 1.0;
    return HPOLYTOPE(d, A, b);
}

static HPOLYTOPE make_cube(int d) {
    MT A = MT::Zero(2 * d, d);
    VT b = VT::Ones(2 * d);
    for (int i = 0; i < d; ++i) {
        A(i, i)     =  1.0; 
        A(d + i, i) = -1.0;   
    }
    return HPOLYTOPE(d, A, b);
}

static MT exact_simplex_samples(int d, int n, unsigned seed) {
    std::mt19937 gen(seed);
    std::exponential_distribution<double> expo(1.0);
    MT S(d, n);
    for (int j = 0; j < n; ++j) {
        std::vector<double> e(d + 1);
        double sum = 0.0;
        for (auto& v : e) { v = expo(gen); sum += v; }
        for (int i = 0; i < d; ++i) S(i, j) = e[i] / sum;
    }
    return S;
}

static MT cube_samples(int d, int n, double last_lo, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> full(-1.0, 1.0);
    std::uniform_real_distribution<double> last(last_lo, 1.0);
    MT S(d, n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < d - 1; ++i) S(i, j) = full(gen);
        S(d - 1, j) = last(gen);
    }
    return S;
}

TEST_CASE("kolmogorov_prob") {
    CHECK(kolmogorov_prob(0.0)  == doctest::Approx(1.0));
    CHECK(kolmogorov_prob(-1.0) == doctest::Approx(1.0));
    CHECK(kolmogorov_prob(0.5)  == doctest::Approx(0.9639452).epsilon(1e-5));
    CHECK(kolmogorov_prob(1.0)  == doctest::Approx(0.2699997).epsilon(1e-5));
    CHECK(kolmogorov_prob(1.36) == doctest::Approx(0.0495).epsilon(2e-2));
    CHECK(kolmogorov_prob(3.0)  == doctest::Approx(0.0).epsilon(1e-6));

    double prev = 1.0;
    for (double z = 0.3; z < 4.0; z += 0.1) {
        double p = kolmogorov_prob(z);
        CHECK(p >= 0.0);
        CHECK(p <= prev + 1e-12);
        prev = p;
    }
}

TEST_CASE("ks_uniform_simplex") {
    const int d = 5, n = 5000;
    HPOLYTOPE P = make_simplex(d);
    MT S = exact_simplex_samples(d, n, 12345);

    auto [ks, p, obs, exp] = global_scaling_test(P, S, 1);

    CHECK(p > 1e-3);         
    CHECK(ks < 0.05);

    REQUIRE(obs.size() == 10);
    REQUIRE(exp.size() == 10);
    for (int k = 0; k < 10; ++k) {
        CHECK(exp[k] == doctest::Approx(0.1 * (k + 1)));
        CHECK(obs[k] == doctest::Approx(exp[k]).epsilon(0.1));
        if (k > 0) CHECK(obs[k] >= obs[k - 1]);            // cumulative
    }
    CHECK(obs[9] == doctest::Approx(1.0));
}

TEST_CASE("ks_nonuniform") {
    const int d = 5, n = 5000;
    HPOLYTOPE P = make_simplex(d);
    MT S = exact_simplex_samples(d, n, 777);

    VT c = VT::Constant(d, 1.0 / (d + 1));
    for (int j = 0; j < n; ++j) S.col(j) = c + 0.6 * (S.col(j) - c);

    auto [ks, p, obs, exp] = global_scaling_test(P, S, 1);
    CHECK(p < 1e-6);
    CHECK(ks > 0.2);
}

TEST_CASE("ks_edge_cases") {
    const int d = 3;
    HPOLYTOPE P = make_simplex(d);
    MT S = exact_simplex_samples(d, 100, 1);

    SUBCASE("too few samples after thinning returns p = 1") {
        auto [ks, p, obs, exp] = global_scaling_test(P, S, 50);
        CHECK(ks == 0.0);
        CHECK(p == 1.0);
        CHECK(obs.empty());
    }
    SUBCASE("thinning_factor < 1 is clamped, not a crash") {
        auto [ks, p, obs, exp] = global_scaling_test(P, S, 0);
        CHECK(p >= 0.0);
        CHECK(p <= 1.0);
    }
}

TEST_CASE("ks_cube_top_half") {
    const int d = 3, n = 5000;
    HPOLYTOPE P = make_cube(d);

    SUBCASE("control: uniform over the whole cube passes") {
        MT S = cube_samples(d, n, -1.0, 42);
        auto [ks, p, obs, expc] = global_scaling_test(P, S, 1);
        CHECK(p > 1e-3);
        CHECK(ks < 0.05);
    }

    SUBCASE("points only in the top half are rejected") {
        MT S = cube_samples(d, n, 0.0, 42);

        for (int j = 0; j < n; ++j) {
            VT x = S.col(j);
            REQUIRE((P.get_mat() * x - P.get_vec()).maxCoeff() <= 1e-12);
        }

        auto [ks, p, obs, expc] = global_scaling_test(P, S, 1);
        CHECK(p < 1e-6);
        CHECK(ks > 0.05);
    }
}
