// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <Eigen/Eigen>
#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/constrained_polytope/bigg_parser.hpp"
#include "preprocess/constrained_polytope/simplifier.hpp"
#include "preprocess/constrained_polytope/transformation.hpp"

typedef double NT;
typedef Cartesian<NT>::Point Point;
typedef ConstrainedPolytope<Point> Polytope;
typedef Polytope::MT MT;
typedef Polytope::VT VT;

NT const INF = std::numeric_limits<NT>::infinity();

// The dimensions every parameterized test runs on
std::vector<unsigned> const DIMS = {2, 5, 10, 20};

// Converts the vector format to a sparse eigen matrix.
MT rows_to_sparse(std::vector<std::vector<NT>> const& rows, unsigned d)
{
    std::vector<Eigen::Triplet<NT>> t;
    for (unsigned i = 0; i < rows.size(); ++i)
        for (unsigned j = 0; j < d; ++j) 
            if (rows[i][j] != NT(0)) t.emplace_back(i, j, rows[i][j]);

    MT A((Eigen::Index)rows.size(), d);
    A.setFromTriplets(t.begin(), t.end());
    return A;
}

VT to_vec(std::vector<NT> const& v) { return Eigen::Map<VT const>(v.data(), v.size()); }

Polytope make_poly(unsigned d,
                   std::vector<std::vector<NT>> const& A_eq, std::vector<NT> const& b_eq,
                   std::vector<std::vector<NT>> const& A_in, std::vector<NT> const& b_in,
                   std::vector<NT> const& lb, std::vector<NT> const& ub)
{
    return Polytope(d, rows_to_sparse(A_eq, d), to_vec(b_eq), rows_to_sparse(A_in, d), to_vec(b_in), to_vec(lb), to_vec(ub));
}

Polytope make_simplex(unsigned d)
{
    std::vector<std::vector<NT>> A_eq(1, std::vector<NT>(d, 1.0));
    return make_poly(d, A_eq, {1.0}, {}, {}, std::vector<NT>(d, 0.0), std::vector<NT>(d, 1.0));
}

// This creates a cube [0,1]^d with redundant rows:
//      x_i <= 2 for every i       (strictly redundant)
//      sum(x) <= d
//      -sum(x) <= 1               (strictly redundant)
//      x_i - x_j <= 1 for i < j
Polytope make_redundant_cube(unsigned d) 
{
    std::vector<std::vector<NT>> A_in;
    std::vector<NT> b_in;

    for (unsigned i = 0; i < d; ++i) {
        std::vector<NT> row(d, 0.0);
        row[i] = 1.0;
        A_in.push_back(row);
        b_in.push_back(2.0);
    }

    A_in.push_back(std::vector<NT>(d, 1.0));
    b_in.push_back(NT(d));

    A_in.push_back(std::vector<NT>(d, -1.0));
    b_in.push_back(1.0);

    for (unsigned i = 0; i < d; ++i) {
        for (unsigned j = i+1; j < d; ++j) {
            std::vector<NT> row(d, 0.0);
            row[i] = 1.0;
            row[j] = -1.0;
            A_in.push_back(row);
            b_in.push_back(1.0);
        }
    }

    return make_poly(d, {}, {}, A_in, b_in, std::vector<NT>(d, 0.0), std::vector<NT>(d, 1.0));
}

// The cube [0,1]^d with the first p = d/2 variables pinned to 1 by a chain of equalities.
// x_0 = 1
// x_j - x_{j-1} = 0 for j = 1, ..., p-1.
// Only x_0 is pinned by a single row, the rest only through the chain.
// The pinned value sits on the upper bound, so the presolve must free the bounds,
// for the interior point to have positive slack.
Polytope make_pinned_cube(unsigned d) 
{
    unsigned const p = d/2;
    std::vector<std::vector<NT>> A_eq;
    std::vector<NT> b_eq;

    std::vector<NT> row(d, 0.0);
    row[0] = 1.0;
    A_eq.push_back(row);
    b_eq.push_back(1.0);

    for (unsigned j = 1; j < p; ++j) {
        std::vector<NT> row(d, 0.0);
        row[j] = 1.0;
        row[j-1] = -1.0;
        A_eq.push_back(row);
        b_eq.push_back(0.0);
    }

    return make_poly(d, A_eq, b_eq, {}, {}, std::vector<NT>(d, 0.0), std::vector<NT>(d, 1.0));
}

// The cube [0,1]^d cut by d-1 facets, each with a parallel redundant copy.
//      x_i + x_{i+1} <= 1.5 for i = 0, ..., d-2   (facet, cuts off the corner x_i = x_{i+1} = 1)
//      x_i + x_{i+1} <= 1.8 for i = 0, ..., d-2   (strictly redundant, implied by the 1.5 row)
// The bounds stay facets: on x_i = 1 the cuts only force x_{i-1}, x_{i+1} <= 0.5.
Polytope make_cut_cube(unsigned d)
{
    std::vector<std::vector<NT>> A_in;
    std::vector<NT> b_in;

    for (NT b : {1.5, 1.8}) {
        for (unsigned i = 0; i+1 < d; ++i) {
            std::vector<NT> row(d, 0.0);
            row[i] = 1.0;
            row[i+1] = 1.0;
            A_in.push_back(row);
            b_in.push_back(b);
        }
    }

    return make_poly(d, {}, {}, A_in, b_in, std::vector<NT>(d, 0.0), std::vector<NT>(d, 1.0));
}

// The cube [0,1]^d with the q = d/2 independent equalities
//      x_{2k} + x_{2k+1} = 1 for k = 0, ..., q-1
// and 2q+1 dependent rows built from them:
//      a copy of every row
//      -10 times every row
//      the sum of all q rows
// The result is a product of q segments and a (d-2q)-cube.
Polytope make_dependent_cube(unsigned d)
{
    unsigned const q = d/2;
    std::vector<std::vector<NT>> A_eq;
    std::vector<NT> b_eq;

    for (unsigned k = 0; k < q; ++k) {
        std::vector<NT> row(d, 0.0);
        row[2*k] = 1.0;
        row[2*k+1] = 1.0;
        A_eq.push_back(row);
        b_eq.push_back(1.0);
    }

    // A copy of every row
    for (unsigned k = 0; k < q; ++k) {
        A_eq.push_back(A_eq[k]);
        b_eq.push_back(b_eq[k]);
    }

    // -10 times every row
    for (unsigned k = 0; k < q; ++k) {
        std::vector<NT> scaled(d, 0.0);
        for (unsigned j = 0; j < d; ++j) scaled[j] = -10*A_eq[k][j];
        A_eq.push_back(scaled);
        b_eq.push_back(-10*b_eq[k]);
    }

    // The sum of all q rows
    std::vector<NT> sum(d, 0.0);
    for (unsigned k = 0; k < q; ++k) 
        for (unsigned j = 0; j < d; ++j)
            sum[j] += A_eq[k][j];

    A_eq.push_back(sum);
    b_eq.push_back(NT(q));

    return make_poly(d, A_eq, b_eq, {}, {}, std::vector<NT>(d, 0.0), std::vector<NT>(d, 1.0));
}

TEST_CASE("simplex") {
    for (unsigned d : DIMS) {
        CAPTURE(d);
        auto res = simplify(make_simplex(d));
        REQUIRE(res.status == SimplifierStatus::OK);

        CHECK(res.polytope.getNumFiniteBounds() == d);

        if (d >= 3) CHECK(res.polytope.getNumFiniteUpperBounds() == 0);

        auto [HP, shift, N] = transform(res.polytope);
        CHECK(HP.dimension() == d-1);
        CHECK(HP.num_of_hyperplanes() == d);
    }
}


TEST_CASE("redundant_cube") {
    for (unsigned d : DIMS) {
        CAPTURE(d);
        auto res = simplify(make_redundant_cube(d));
        REQUIRE(res.status == SimplifierStatus::OK);

        CHECK(res.polytope.getNumInequalities() == 0);
        CHECK(res.polytope.getNumFiniteLowerBounds() == d);
        CHECK(res.polytope.getNumFiniteUpperBounds() == d);

        auto [HP, shift, N] = transform(res.polytope);
        CHECK(HP.dimension() == d);
        CHECK(HP.num_of_hyperplanes() == 2*d);
    }
}

TEST_CASE("pinned_variables") {
    for (unsigned d : DIMS) {
        CAPTURE(d);
        unsigned const p = d/2;
        auto res = simplify(make_pinned_cube(d));
        REQUIRE(res.status == SimplifierStatus::OK);

        CHECK(res.presolve.pinned_vars == p);
        CHECK(res.presolve.rank == p);

        // The pinned vars need to lose their bounds, the free ones should not
        auto const& lb = res.polytope.getLowerBounds();
        auto const& ub = res.polytope.getUpperBounds();

        for (unsigned j = 0; j < d; ++j) {
            CAPTURE(j);
            CHECK(std::isfinite(lb(j)) == (j >= p));
            CHECK(std::isfinite(ub(j)) == (j >= p));
        }

        // Check that the interior point has the pinned values.
        REQUIRE(res.interior_point.valid);
        for (unsigned j = 0; j < p; ++j)
            CHECK(res.interior_point.point(j) == doctest::Approx(1.0));

        // The affine hull should be a d-p dimensional cube
        auto [HP, shift, N] = transform(res.polytope);
        CHECK(HP.dimension() == d-p);
        CHECK(HP.num_of_hyperplanes() == 2*(d-p));
    }
}

TEST_CASE("keep_facets") {
    for (unsigned d : DIMS) {
        CAPTURE(d);
        auto res = simplify(make_cut_cube(d));
        REQUIRE(res.status == SimplifierStatus::OK);

        // The d-1 cuts stay, the copies go
        CHECK(res.polytope.getNumInequalities() == d-1);
        auto const& b_in = res.polytope.getInequalityRHS();
        for (unsigned i = 0; i < res.polytope.getNumInequalities(); ++i)
            CHECK(b_in(i) == doctest::Approx(1.5));

        // Every bound stays
        CHECK(res.polytope.getNumFiniteLowerBounds() == d);
        CHECK(res.polytope.getNumFiniteUpperBounds() == d);

        auto [HP, shift, N] = transform(res.polytope);
        CHECK(HP.dimension() == d);
        CHECK(HP.num_of_hyperplanes() == 3*d-1);
    }
}

TEST_CASE("dependent_equalities") {
    for (unsigned d : DIMS) {
        CAPTURE(d);
        unsigned const q = d/2;
        auto res = simplify(make_dependent_cube(d));
        REQUIRE(res.status == SimplifierStatus::OK);

        // The input has 3q+1 equalities of rank q
        CHECK(make_dependent_cube(d).getNumEqualities() == 3*q+1);

        // The 2q+1 dependent rows go, the q independent ones stay
        CHECK(res.presolve.rank == q);
        CHECK(res.dependent_rows.rank == q);
        CHECK(res.dependent_rows.removed_rows == 2*q+1);
        CHECK(res.polytope.getNumEqualities() == q);

        // No variable should be pinned.
        CHECK(res.presolve.pinned_vars == 0);
        CHECK(res.reduce.fixed_vars == 0);

        // q segments and a (d-2q)-dimensional cube
        auto [HP, shift, N] = transform(res.polytope);
        CHECK(HP.dimension() == d-q);
        CHECK(HP.num_of_hyperplanes() == 2*q+2*(d-2*q));

    }
}

TEST_CASE("zero_reduction_passes") {
    SimplifierConfig c;
    c.max_reduction_passes = 0;

    // Has a strict interior point, so Clarkson runs without any reduction
    auto ok = simplify(make_simplex(3), c);
    CHECK(ok.status == SimplifierStatus::OK);
    CHECK(ok.reduce.passes == 0);

    // x_2 is fixed by lb = ub, so only the reduce stage can find it
    auto bad = simplify(make_poly(3, {}, {}, {}, {}, {0, 0, 2}, {1, 1, 2}), c);
    CHECK(bad.status == SimplifierStatus::NO_INTERIOR_POINT);
}

