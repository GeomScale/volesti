// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>

#include "random.hpp"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"
#include "optimization/simulated_annealing.hpp"


typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;


void call_test_solve_sdp() {

    // read the sdp
    std::string fileName("sdp__2_8.txt");
    std::cout << fileName <<"\n\n";fflush(stdout);
    SPECTRAHEDRON spectrahedron;
    Point objFunction;

    std::ifstream in;
    in.open(fileName, std::ifstream::in);
    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.loadSDPAFormatFile(in, spectrahedron, objFunction);

    // We will need an initial interior point. In this
    // spectrahedron the origin (zero point) is interior
    Point initialPoint(spectrahedron.getLMI().dimension());

    // First some parameters for the solver
    // desired relative error
    NT rel_error = 0.001;

    // Declare settings
    SimulatedAnnealingSettings<Point> settings(rel_error);



    bool correct = false;
    int tries = 0;

    // try some times; since it is randomized it may fail...
    while (!correct && tries < 5){
        // solve the program
        Point sol;
        NT min = solve_sdp(spectrahedron, objFunction, settings, initialPoint, sol);
        // assume optimal solution = -1.3888
        double relativeError = std::fabs((-1.3888 - min) / -1.3888);
        correct = relativeError <= 0.01;
        tries++;
    }

    CHECK(correct);
}

TEST_CASE ("test_solve_sdp") {
    call_test_solve_sdp(); // test the simulated annealing SDP solver
}