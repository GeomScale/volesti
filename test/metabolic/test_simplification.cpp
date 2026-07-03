// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include <limits>
#include <vector>
#include "doctest.h"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/metabolic/simplification_exhaustive.hpp"
#include "preprocess/metabolic/simplification_clarkson.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;
typedef typename Polytope::MT MT;
typedef typename Polytope::VT VT;
static NT const INF = std::numeric_limits<NT>::infinity();

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

unsigned bounds_relaxed(Polytope const& Po, Polytope const& Ps) {
    return Po.getNumFiniteBounds()-Ps.getNumFiniteBounds();
}

unsigned dims_fixed(Polytope const& Po, Polytope const& Ps) {
    return Ps.getNumEqualities()-Po.getNumEqualities();
}

template <typename Simplifier>
void test_cube_no_change(unsigned d) 
{   
    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope::cube(d);
    auto [P2, success] = Simplifier::run(P1, false);

    REQUIRE(success);
    CHECK(bounds_relaxed(P1, P2) == 0);
    CHECK(dims_fixed(P1, P2) == 0);
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
        b_l(i) = -10.0;
        b_u(i) = 10.0;
    }

    for (unsigned i = d; i < m; ++i) {
        b_l(i) = 0.0;
        b_u(i) = 1.0;
    }

    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope(m, A_eq, b_l, b_u, b_eq);
    auto [P2, success] = Simplifier::run(P1, false);

    REQUIRE(success);

    // Checks the simplification statistics.
    CHECK(bounds_relaxed(P1, P2) == 2*d);
    CHECK(dims_fixed(P1, P2) == 0);

    // Checks that equalities were untouched.
    CHECK(P2.getEqualities().isApprox(P1.getEqualities()));
    CHECK(P2.getEqualityBounds() == P1.getEqualityBounds());

    // Checks that the first d reactions were relaxed.
    CHECK((P2.getLowerBounds().head(d).array() == -INF).all());
    CHECK((P2.getUpperBounds().head(d).array() == INF).all());

    // The other m-d bounds must remain untouched.
    CHECK(P2.getLowerBounds().tail(d) == P1.getLowerBounds().tail(d));
    CHECK(P2.getUpperBounds().tail(d) == P1.getUpperBounds().tail(d));
}

template <typename Simplifier>
void test_simplex_relaxed_bounds(unsigned d) 
{
    INFO("simplifier: " << Simplifier::name());

    Polytope P1 = Polytope::simplex(d);
    auto [P2, success] = Simplifier::run(P1, false);

    REQUIRE(success);

    // Checks statistics.
    CHECK(bounds_relaxed(P1, P2) == d);
    CHECK(dims_fixed(P1, P2) == 0);

    // Checks every upper bound is infinity.
    CHECK((P2.getUpperBounds().array() == INF).all());

    // Checks that the rest are untouched.
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
    auto [P2, success] = Simplifier::run(P1, true);

    REQUIRE(success);

    // Checks statistics.
    CHECK(bounds_relaxed(P1, P2) == 2*d);
    CHECK(dims_fixed(P1, P2) == d);

    // Checks that the pinned constraints were added to A_eq.
    CHECK(P2.getNumEqualities() == d);
    CHECK((unsigned)P2.getEqualityBounds().size() == d);

    // Checks that the last d reactions were relaxed.
    CHECK((P2.getLowerBounds().tail(d).array() == -INF).all());
    CHECK((P2.getUpperBounds().tail(d).array() == INF).all());

    // Checks that the rest of the bounds remain untouched.
    CHECK(P2.getLowerBounds().head(d) == P1.getLowerBounds().head(d));
    CHECK(P2.getUpperBounds().head(d) == P1.getUpperBounds().head(d));
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
