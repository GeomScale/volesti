// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "Eigen/Eigen"
#include "vector"
#include <fstream>
#include "random.hpp"
#include "cartesian_geom/cartesian_kernel.h"
#include "spectrahedron.h"
#include "SDPAFormatManager.h"
#include "string"
#include "iostream"

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<NT, MT, VT> SPECTRAHEDRON;

// create the data to write to sdpa file
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


void call_test_sdpa_format() {

    //create data to write to file
    SPECTRAHEDRON spectrahedron;
    Point objFunction;

    populateValues(spectrahedron, objFunction);

    std::string fileName("sdpaFormatTest.txt");

    // write the data to a file
    std::filebuf fb;
    fb.open(fileName, std::ios::out);
    std::ostream os(&fb);
    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.writeSDPAFormatFile(os, spectrahedron, objFunction);
    os.flush();

    // now Read it back!
    SPECTRAHEDRON spectrahedron1;
    Point objFunction1;
    std::ifstream in;
    in.open(fileName, std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in, spectrahedron1, objFunction1);

    // check if what you wrote and read are the same
    bool areTheSame = true;

    while (true) {
        if (!(objFunction == objFunction1)) {
            areTheSame = false;
            break;
        }

        std::vector<MT> matrices1 = spectrahedron1.getLMI().getMatrices();

        int atMatrix1 = 0;
        std::vector<MT> const & _matrices = spectrahedron.getLMI().getMatrices();

        for (auto matrix = _matrices.begin(); matrix != _matrices.end() ; matrix++) {
            MT matrix1 = matrices1[atMatrix1++];

            for (int i = 0; i < matrix->rows(); ++i) {
                for (int j = 0; j < matrix->cols(); ++j) {
                    if ((*matrix)(i, j) != matrix1(i, j)) {
                        areTheSame = false;
                        break;
                    }
                }
            }
        }

        break;
    }


    CHECK(areTheSame);
}


TEST_CASE ("sdpa_format_parser") {
    call_test_sdpa_format();
}

