// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file


// Edited by HZ on 11.06.2020 - mute doctest.h
#include "doctest.h"
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
NT factorial(NT n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT>
void test_values(NT volume, NT expected, NT exact)
{
    std::cout << "Computed volume " << volume << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    std::cout << "Relative error (expected) = "
              << std::abs((volume-expected)/expected) << std::endl;
    std::cout << "Relative error (exact) = "
              << std::abs((volume-exact)/exact) << std::endl;
    CHECK((std::abs((volume - exact)/exact) < 0.35 || 
           std::abs((volume - expected)/expected) < 0.00001));
}

/*
template <class Polytope>
void test_volume(Polytope &HP,
                 double const& expectedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 10 + HP.dimension()/10;
    NT e=0.1, volume;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    volume = volume_cooling_balls<BilliardWalk, RNGType>(HP, , walk_len, e).second;
    test_values(volume, expectedBilliard, exact);
}*/

template <typename NT>
void call_test_volume(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef Spectrahedron<Point> spectrahedron;

    SdpaFormatManager<NT> sdpaFormatManager;

    std::ifstream in1;
    spectrahedron spectra;
    Point objFunction;
    in1.open("spectra_data/sdp_prob_20_20.txt", std::ifstream::in);
    sdpaFormatManager.loadSDPAFormatFile(in1, spectra, objFunction);

    // Setup the parameters
    int walk_len = 10 + spectra.dimension()/10;
    NT e = 0.1, volume;
    Point interior_point(spectra.dimension());
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    
    std::cout << "--- Testing spectrahedron 20_20" << std::endl;
    volume = volume_cooling_balls<BilliardWalk, RNGType>(spectra, interior_point, walk_len, e).second;
    //test_volume(spectra, 1118.63, 1118.63);
    test_values(volume, 1118.63, 1118.63);
}


TEST_CASE("spectra_volume") {
    call_test_volume<double>();
}


