// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include <vector>
#include "doctest.h"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/metabolic_polytope.hpp"
#include "preprocess/metabolic/scaling.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;
typedef typename Polytope::MT MT;
typedef typename Polytope::VT VT;

// Builds the polytope:
//
// 1e-5*v0 - v2 = 0
// 1000*v1 - v3 = 0
Polytope build_polytope() {
    unsigned d = 4, m = 2;
    NT const INF = std::numeric_limits<NT>::infinity();

    std::vector<Polytope::Triplet> triplets;
    triplets.emplace_back(0, 0, NT(1e-5));
    triplets.emplace_back(0, 2, NT(-1));
    triplets.emplace_back(1, 1, NT(1000));
    triplets.emplace_back(1, 3, NT(-1));

    MT A_eq(m, d);
    A_eq.setFromTriplets(triplets.begin(), triplets.end());
    A_eq.makeCompressed();

    VT b_eq = VT::Zero(m);

    VT b_l(d), b_u(d);
    b_l << -1e3, -1e-3, -INF, -1.0;
    b_u <<  1e3,  1e-3, INF, 1.0;

    return Polytope(d, A_eq, b_l, b_u, b_eq);
}

// Builds points of the polytope above.
std::vector<VT> build_points() {
    std::vector<VT> points;
    double v0[] = {0.0, 1.0, -1.0, 100.0, -100.0, 1e3, -1e3};
    double v1[] = {0.0, 1e-6, -1e-6, 1e-4, -1e-4, 1e-3, -1e-3};

    for (double a : v0) {
        for (double b : v1) {
            VT x(4);
            x << a, b, 1e-5*a, 1000.0*b;
            points.push_back(x);
        }
    }

    return points;
}

void test_no_scaling(Polytope const& P) {
    Scaling<Point> s;
    Polytope SP = scale(P, s, NoScaling{});
    CHECK(s.col.isOnes());
    CHECK(s.row.isOnes());
}

template <typename ScalingPolicy>
void test_scaling(Polytope const& P, std::vector<VT> const& points) {
    Scaling<Point> s;
    Polytope SP = scale(P, s, ScalingPolicy{});

    // Scaling dimensions must match.
    REQUIRE((unsigned)s.col.size() == P.getDimension());
    REQUIRE((unsigned)s.row.size() == (unsigned)P.getEqualities().rows());

    // Factors must be positive/finite.
    CHECK((s.col.array() > 0.0).all());
    CHECK((s.row.array() > 0.0).all());
    CHECK(s.col.allFinite());
    CHECK(s.row.allFinite());

    // Dimensions after must match.
    CHECK(SP.getDimension() == P.getDimension());
    CHECK(SP.getNumFiniteBounds() == P.getNumFiniteBounds());

    // Scales back the polytope, this should return P.
    Polytope OP = rescale(SP, s, true);

    CHECK(OP.getLowerBounds() == P.getLowerBounds());
    CHECK(OP.getUpperBounds() == P.getUpperBounds());
    CHECK(OP.getEqualityBounds() == P.getEqualityBounds());
    CHECK((OP.getEqualities()-P.getEqualities()).norm() == doctest::Approx(0.0));

    for (VT const& x : points) {
        VT y = scale_point(x, s, true);
        VT residual = SP.getEqualities()*y-SP.getEqualityBounds();

        // Checks that the scaled point belong to the scaled polytope.
        CHECK(residual.cwiseAbs().maxCoeff() == doctest::Approx(0.0));
        CHECK((y.array() >= SP.getLowerBounds().array()).all());
        CHECK((y.array() <= SP.getUpperBounds().array()).all());

        // Rescales the point and checks that it is the original point.
        VT z = scale_point(y, s, false);
        CHECK(x == z);
    }
}

TEST_CASE("test_no_scaling"){
    test_no_scaling(build_polytope());
}

TEST_CASE_TEMPLATE("test_scaling", ScalingPolicy, NoScaling, MaxBoundScaling, GMScaling) {
    test_scaling<ScalingPolicy>(build_polytope(), build_points());
}