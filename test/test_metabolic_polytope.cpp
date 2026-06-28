// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/metabolic_polytope.h"
#include "io/bigg_parser.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;
typedef typename Polytope::MT MT;
typedef typename Polytope::VT VT;

static const std::string ECOLI_JSON = std::string(BIGG_DIR)+"/e_coli_core.json";

void test_construction(unsigned d) {
    VT b_u = VT::Ones(d);
    VT b_l = -VT::Ones(d);
    MT A_eq(0, d);
    VT b_eq = VT::Zero(0);

    Polytope P(d, A_eq, b_l, b_u, b_eq);

    CHECK(P.getDimension() == d);
    CHECK(P.getEqualityBounds() == b_eq);
    CHECK(P.getEqualities().isApprox(A_eq));
    CHECK(P.getLowerBounds() == b_l);
    CHECK(P.getUpperBounds() == b_u);
    CHECK(P.getNumEqualities() == 0);
    CHECK(P.getNumFiniteBounds() == 2*d);
}

void test_copy(unsigned d) {
    Polytope P1 = Polytope::cube(d);
    Polytope P2(P1);

    CHECK(P1.getDimension() == P2.getDimension());
    CHECK(P1.getEqualityBounds() == P2.getEqualityBounds());
    CHECK(P1.getEqualities().isApprox(P2.getEqualities()));
    CHECK(P1.getLowerBounds() == P2.getLowerBounds());
    CHECK(P1.getUpperBounds() == P2.getUpperBounds());
}

void test_ecoli_construction() {
    Polytope P = bigg::parse_from_json<Point>(ECOLI_JSON);

    CHECK(P.getDimension() == 95);
    CHECK(P.getEqualities().rows() == 72);
    CHECK(P.getEqualities().cols() == 95);
    CHECK(P.getLowerBounds().size() == 95);
    CHECK(P.getUpperBounds().size() == 95);
    CHECK(P.getEqualityBounds().size() == 72);
    CHECK(P.getEqualityBounds().isZero());
    CHECK((P.getLowerBounds().array() <= P.getUpperBounds().array()).all());
}

TEST_CASE("test_construction") {
    test_construction(10);
    test_ecoli_construction();
    test_copy(10);
}