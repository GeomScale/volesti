// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022 Alexandros Manochis

// Licensed under GNU LGPL.3, see LICENCE file

#include <fstream>
#include <iostream>
#include "misc.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

#include "matrix_operations/EigenvaluesProblems.h"
#include "SDPAFormatManager.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"


template <typename NT>
void example_volume(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef Spectrahedron<Point> spectrahedron;

    SdpaFormatManager<NT> sdpaFormatManager;

    std::ifstream in1;
    spectrahedron spectra;
    Point objFunction;
    in1.open("./../../test/spectra_data/sdp_prob_20_20.txt", std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in1, spectra, objFunction);

    // Setup the parameters
    int walk_len = 2;
    NT e = 0.1, volume;
    Point interior_point(spectra.dimension());
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    
    std::cout << "--- Testing spectrahedron 20_20" << std::endl;
    volume = volume_cooling_balls<AcceleratedBilliardWalk, RNGType>(spectra, interior_point, walk_len, e).second;
    
    std::cout<<"Computed volume = "<<volume<<std::endl;
}


int main() {

    example_volume<double>();

    return 0;
}


