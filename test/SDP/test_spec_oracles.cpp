// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "Eigen/Eigen"
#include <vector>
#include "random.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "spectrahedron.h"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<NT, MT, VT> SPECTRAHEDRON;


// create a spectrahedron to test the boundary oracles
void populateValues(SPECTRAHEDRON &spectrahedron, Point &objFunction) {
    VT c(2);
    c(0) = 1;
    c(1) = 1;
    objFunction = Point(c);

    MT A0, A1, A2;
    A0.setZero(3, 3);
    A1.setZero(3, 3);
    A2.setZero(3, 3);

    A0(0, 0) = -1;
    A0(1, 1) = -2;
    A0(1, 2) = A0(2, 1) = 1;
    A0(2, 2) = -2;

    A1(0, 0) = -1;
    A1(1, 2) = A1(2, 1) = 1;

    A2(0, 2) = A2(2, 0) = -1;

    std::vector<MT> matrices;
    matrices.push_back(A0);
    matrices.push_back(A1);
    matrices.push_back(A2);

    LMI<NT, MT, VT> lmi(matrices);
    spectrahedron = SPECTRAHEDRON(lmi);
}


void call_test_coordinate_intersection() {

    // create a spectrahedron to test the boundary oracles
    SPECTRAHEDRON spectrahedron;
    Point objFunction;
    populateValues(spectrahedron, objFunction);

    // create an initial point in the spectrahedron
    VT initialPoint;
    initialPoint.setZero(2);

    SPECTRAHEDRON::PrecomputedValues precomputedValues;
    std::pair<NT, NT> ret = spectrahedron.coordinateIntersection(initialPoint, 2, precomputedValues);

    bool correct = true;

    // actual result is
    // 1.22474 -1.22474

    if (std::fabs(ret.first - 1.22474) > 0.0001 || std::fabs(ret.second + 1.22474) > 0.0001)
        correct = false;

    CHECK(correct);
}

void call_test_positive_intersection() {

    // we will test with the curve at^2 + bt + c
    // where c = interior point
    // and a = objective function

    SPECTRAHEDRON spectrahedron;
    Point a;
    populateValues(spectrahedron, a);

    // create an initial point in the spectrahedron
    // the origin will do!
    VT c;
    c.setZero(2);

    // create a b
    VT b;
    b.setZero(2);
    b(0) = 1;
    b(1) = 1;

    SPECTRAHEDRON::PrecomputedValues precomputedValues;
    NT ret = spectrahedron.positiveIntersection(a.getCoefficients(), b, c, precomputedValues);

    bool correct = true;

    // actual result is
    // ret = 0.5294590425

    if (std::fabs(ret - 0.5294590425) > 0.0001)
        correct = false;

    CHECK(correct);
}


TEST_CASE ("test_spec_oracles") {
    call_test_positive_intersection(); // hmc boundary oracle
    call_test_coordinate_intersection(); // coordinate hit and run oracle
}

