// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

// This examples illustrates how to sample a spectrahedron under the Boltzmann distribution with
// HMC random walk. It will read the spectrahedron from data/sdp_n2m3.txt.

//#define VOLESTI_DEBUG


#include <iostream>
#include <fstream>

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"
#include "random_walks/boltzmann_hmc_walk.hpp"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef BoltzmannHMCWalk::Walk<SPECTRAHEDRON, RNGType> HmcWalk;
typedef BoltzmannHMCWalk::Walk<SPECTRAHEDRON, RNGType>::Settings HmcWalkSettings;


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

    // required parameters for the random walk
    int walkLength = 5;
    RNGType randomNumberGenerator(spectrahedron.getLMI().dimension()); // this class provides random numbers
    NT temperature = 1;

    // estimate the diameter of the body
    int pointsNum = 10;
    NT diameter = spectrahedron.estimateDiameter(pointsNum, initialPoint);

    // declare the settings and
    HmcWalkSettings settings(walkLength, randomNumberGenerator, objFunction, temperature, diameter);

    // declare the random walk
    HmcWalk hmcWalk(settings);

    // sample three points from the spectrahedron
    pointsNum = 3;
    std::list<Point> points;
    hmcWalk.apply(spectrahedron, initialPoint, pointsNum, points);

    // print sampled points
    for (Point point : points)
        point.print();

    return 0;
}

