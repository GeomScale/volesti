// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/metabolic/exhaustive_simplification.hpp"
#include "preprocess/metabolic/clarkson_simplification.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;
typedef typename Polytope::MT MT;
typedef typename Polytope::VT VT;

struct Exhaustive {
    static exhaustive_simplification::Result<Point> run(Polytope const& P, 
                                                        bool fix_dimensions)
    {
        exhaustive_simplification::Config config;
        config.fix_dimensions = fix_dimensions;
        return exhaustive_simplification::simplify(P, config);
    }

    static char const* name() {return "exhaustive";}
};

struct Clarkson {
    static clarkson_simplification::Result<Point> run(Polytope const& P, 
                                                      bool fix_dimensions)
    {
        clarkson_simplification::Config config;
        config.fix_dimensions = fix_dimensions;
        return clarkson_simplification::simplify(P, config);
    }

    static char const* name() {return "clarkson";}
};

template <typename Simplifier>
void test_cube_no_change(unsigned d) 
{   
    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope::cube(d);
    auto result = Simplifier::run(P1, false);
    Polytope P2 = result.P;

    CHECK(result.bounds_relaxed == 0);
    CHECK(result.dims_fixed == 0);
    CHECK(P2.getEqualities().isApprox(P1.getEqualities()));
    CHECK(P2.getLowerBounds() == P1.getLowerBounds());
    CHECK(P2.getUpperBounds() == P1.getUpperBounds());
    CHECK(P2.getEqualityBounds() == P1.getEqualityBounds());
}

template <typename Simplifier>
void test_cube_relaxed_bounds(unsigned d) 
{
    unsigned m = 2*d;
    VT b_l(m);
    VT b_u(m);
    MT A_eq(d, m);
    VT b_eq = VT::Ones(d);

    std::vector<Polytope::Triplet> triplets;
    for (unsigned i = 0; i < d; ++i) {
        triplets.emplace_back(i, i, NT(1));
        triplets.emplace_back(i, d+i, NT(1));
    }
    A_eq.setFromTriplets(triplets.begin(), triplets.end());
    A_eq.makeCompressed();

    for (unsigned i = 0; i < d; ++i) {
        b_l(i) = - 10.0;
        b_u(i) = 10.0;
    }

    for (unsigned i = d; i < m; ++i) {
        b_l(i) = 0.0;
        b_u(i) = 1.0;
    }

    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope(m, A_eq, b_l, b_u, b_eq);
    auto result = Simplifier::run(P1, false);
    Polytope P2 = result.P;

    CHECK(result.bounds_relaxed == 2*d);
    CHECK(result.dims_fixed == 0);
    CHECK(P2.getEqualities().isApprox(P1.getEqualities()));
    CHECK(P2.getEqualityBounds() == P1.getEqualityBounds());
}

template <typename Simplifier>
void test_simplex_relaxed_bounds(unsigned d) 
{
    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope::simplex(d);
    auto result = Simplifier::run(P1, false);
    Polytope P2 = result.P;

    // Every upper bound is infinity
    for (unsigned i = 0; i < d; ++i)
        CHECK(std::isinf((double)P2.getUpperBounds()(i)));

    CHECK(result.bounds_relaxed == d);
    CHECK(result.dims_fixed == 0);
    CHECK(P1.getDimension() == P2.getDimension());
    CHECK(P1.getLowerBounds() == P2.getLowerBounds());
    CHECK(P1.getEqualities().isApprox(P2.getEqualities()));
    CHECK(P1.getEqualityBounds() == P2.getEqualityBounds());
}

template <typename Simplifier>
void test_cube_degenerate_dimensions(unsigned d) 
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
    auto result = Simplifier::run(P1, true);

    CHECK(result.bounds_relaxed == 2*d);
    CHECK(result.dims_fixed == d);
}

TEST_CASE_TEMPLATE("test_no_change", Simplifier, Exhaustive, Clarkson) {
    test_cube_no_change<Simplifier>(10);
}

TEST_CASE_TEMPLATE("test_no_dimension_fixing", Simplifier, Exhaustive, Clarkson) {
    test_simplex_relaxed_bounds<Simplifier>(10);
    test_cube_relaxed_bounds<Simplifier>(10);
}

TEST_CASE_TEMPLATE("test_dimension_fixing", Simplifier, Exhaustive, Clarkson) {
    test_cube_degenerate_dimensions<Simplifier>(10);
}
