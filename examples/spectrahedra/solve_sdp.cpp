// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

// This examples illustrates how to solve a semidefinite program.
// It will read a semidefinite program from data/sdp_n2m3.txt, solve it and print its solution (minimal value).


//#define VOLESTI_DEBUG
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
typedef Spectrahedron <Point> SPECTRAHEDRON;


int main(int argc, char* argv[]) {
    std::string fileName("data/sdp_n2m3.txt");
    std::string outputFile("new_sdp_n2m3.txt");

    SPECTRAHEDRON spectrahedron;
    Point objFunction;

    // read the spectrahedron
    // open a stream to read the input file
    std::ifstream in;
    in.open(fileName, std::ifstream::in);

    // read the file
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

    // solve the program
    Point sol;
    bool verbose = true;
    NT min = solve_sdp(spectrahedron, objFunction, settings, initialPoint, sol ,verbose);

    // print solution
    std::cout << min << "\n" << "point: ";
    sol.print();

    return 0;
}

