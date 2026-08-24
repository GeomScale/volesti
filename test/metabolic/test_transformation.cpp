// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include "doctest.h"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/metabolic/simplification_exhaustive.hpp"
#include "preprocess/metabolic/simplification_clarkson.hpp"
#include "preprocess/metabolic/transformation.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/boost_random_number_generator.hpp"
#include "sampling/sampling.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;
typedef HPolytope<Point> Hpolytope;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DenseMT;
typedef typename Polytope::MT MT;
typedef typename Polytope::VT VT;
typedef BoostRandomNumberGenerator<boost::mt19937, double> RNG;

struct Exhaustive {
    static std::pair<Polytope, bool> run(Polytope const& P, 
                                         bool fix_dimensions)
    {
        ExhaustiveConfig config;
        config.fix_dimensions = fix_dimensions;
        ExhaustiveSimplifier<Point> simplifier(P, config);
        return simplifier.simplify();
    }

    static std::string name() {return "exhaustive";}
};

struct Clarkson {
    static std::pair<Polytope, bool> run(Polytope const& P, 
                                         bool fix_dimensions)
    {
        ClarksonConfig config;
        config.fix_dimensions = fix_dimensions;
        ClarksonSimplifier<Point> simplifier(P, config);
        return simplifier.simplify();
    }

    static std::string name() {return "clarkson";}
};

NT compute_median_volume(HPolytope<Point> & HP,
                         double e,
                         unsigned walk_len,
                         unsigned num_trials) 
{
    std::vector<double> volumes;
    for (unsigned i = 0; i < num_trials; ++i) {
        RNG rng(HP.dimension());
        rng.set_seed(i);
        auto v = volume_cooling_balls<BallWalk, HPolytope<Point>>(HP, rng, e, walk_len);
        volumes.push_back(v.second);
    }             

    std::stable_sort(volumes.begin(), volumes.end());
    return volumes[volumes.size()/2];
}

bool is_feasible(Polytope const& P, VT const& x, double tol = 1e-10) {
    MT const& A_eq = P.getEqualities();
    VT const& b_eq = P.getEqualityBounds();
    VT const& b_l = P.getLowerBounds();
    VT const& b_u = P.getUpperBounds();

    // Checks that it meets lower/upper bounds
    for (unsigned i = 0; i < (unsigned)x.size(); ++i) {
        if (!std::isinf((double)b_l(i)) && x(i) < b_l(i)-tol) return false;
        if (!std::isinf((double)b_u(i)) && x(i) > b_u(i)+tol) return false;
    }

    // Checks that it meets equalities
    if (!A_eq.rows()) return true;
    VT val = A_eq*x-b_eq;
    return val.cwiseAbs().maxCoeff() < tol; // basically checks that x solves all equations
}

bool check_sampling(Polytope const& P,
                    Hpolytope & HP,
                    DenseMT const& N,
                    VT const& shift,
                    unsigned samples) 
{
    Point c = HP.ComputeInnerBall().first;
    RNG rng(HP.dimension());
    rng.set_seed(0);
    std::list<Point> sampled;
    uniform_sampling<BilliardWalk>(
        sampled, HP, rng, 1, samples, c, 0
    );

    for (auto const& pt : sampled) {
        VT x = N*pt.getCoefficients()+shift;
        if (!is_feasible(P, x)) return false;
    }
    return true;
}

template <typename Simplifier>
void test_cube_transformation(unsigned d) 
{
    unsigned m = 2*d;
    VT b_l(m);
    VT b_u(m);
    MT A_eq(0, m);
    VT b_eq(0);

    for (unsigned i = 0; i < d; ++i) {
        b_l(i) = 0.0;
        b_u(i) = 1.0;
    }

    for (unsigned i = d; i < m; ++i) {
        b_l(i) = 1.0;
        b_u(i) = 1.0+1e-12;
    }

    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope(m, A_eq, b_l, b_u, b_eq);
    auto [P2, ok] = Simplifier::run(P1, true);
    
    REQUIRE(ok);

    auto trs_res = transform(P2);
    HPolytope<Point> HP = std::get<0>(trs_res);
    VT shift = std::get<1>(trs_res);
    DenseMT N = std::get<2>(trs_res);
    unsigned walk_len = 10+d/10;
    double volume = compute_median_volume(HP, 0.05, walk_len, 10);

    CHECK(HP.dimension() == d);
    CHECK(volume > 0.90);                        // checks volume
    CHECK(volume < 1.10);
    CHECK(check_sampling(P1, HP, N, shift, 25)); // checks sampling
}

template <typename Simplifier>
void test_simplex_transformation(unsigned d) 
{   
    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope::simplex(d);
    auto [P2, ok] = Simplifier::run(P1, true);
    
    REQUIRE(ok);

    auto trs_res = transform(P2);
    HPolytope<Point> HP = std::get<0>(trs_res);
    VT shift = std::get<1>(trs_res);
    DenseMT N = std::get<2>(trs_res);

    unsigned walk_len = 10+d/10;
    double actual_v = std::sqrt(d)/std::tgamma(d);
    double volume = compute_median_volume(HP, 0.05, walk_len, 10);

    CHECK(HP.dimension() == d-1);
    CHECK(volume > 0.90*actual_v);               // checks volume
    CHECK(volume < 1.10*actual_v);
    CHECK(check_sampling(P1, HP, N, shift, 25)); // checks sampling
}

TEST_CASE_TEMPLATE("test_cube_transformation", Simplifier, Exhaustive, Clarkson) {
    test_cube_transformation<Simplifier>(10);
}

TEST_CASE_TEMPLATE("test_simplex_transformation", Simplifier, Exhaustive, Clarkson) {
    test_simplex_transformation<Simplifier>(10);
}