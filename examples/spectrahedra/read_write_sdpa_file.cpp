// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

// This examples illustrates how to read and write SDPA format files.
// It will read a semidefinite program from data/sdp_n2m3.txt, print it and then write it to a new file

#define VOLESTI_DEBUG

#include <fstream>
#include <iostream>

#include "random.hpp"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;


int main(int argc, char* argv[]) {
    std::string fileName("data/sdp_n2m3.txt");
    std::string outputFile("new_sdp_n2m3.txt");

    SPECTRAHEDRON spectrahedron;
    Point objFunction;

    // read the semidefinite program
    // and create a vector (objective function) and a spectrahedron

    // open a stream to read the input file
    std::ifstream in;
    in.open(fileName, std::ifstream::in);

    // read the file
    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.loadSDPAFormatFile(in, spectrahedron, objFunction);

    // print the contents
    std::cout << "The objective Function:\n\n";
    objFunction.print();
    std::cout << "\n\nThe matrices of the spectrahedron:\n\n";
    spectrahedron.getLMI().print();

    // open a stream to an output file
    std::filebuf fb;
    fb.open(outputFile, std::ios::out);
    std::ostream os(&fb);

    // write a SDPA format file using the data we read before
    sdpaFormatManager.writeSDPAFormatFile(os, spectrahedron, objFunction);

    return 0;
}

