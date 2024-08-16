// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis
// Copyright (c) 2024 Vladimir Necula

// Contributed and/or modified by Vladimir Necula, as part of Google Summer of
// Code 2024 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_nonspherical_gaussians_crhmc.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include "misc/misc.h"

const unsigned int FIXED_SEED = 42;

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt11213b, NT, FIXED_SEED> RandomNumberGenerator;
typedef boost::mt19937 PolyRNGType;
typedef HPolytope<Point> HPOLYTOPE;

NT nonspherical_crhmc_volume(HPOLYTOPE& polytope) {
    int walk_len = 10;
    NT e = 0.1;
    RandomNumberGenerator rng;
    // NT volume = volume_cooling_gaussians<GaussianBallWalk, RandomNumberGenerator>(polytope, e, walk_len);
    NT volume = non_spherical_crhmc_volume_cooling_gaussians<HPOLYTOPE, RandomNumberGenerator>(polytope, rng, e, walk_len);
    return volume;
}

NT spherical_gaussians_volume(HPOLYTOPE& polytope) {
    int walk_len = 10;
    NT e = 0.1;
    RandomNumberGenerator rng;
    NT volume = volume_cooling_gaussians<GaussianCDHRWalk, RandomNumberGenerator>(polytope, e, walk_len);
    return volume;
}

int main() {

    HPOLYTOPE cube3 = generate_cube<HPOLYTOPE>(3, false);
    std::cout << "Cube3 \n";
    std::cout << "Calculated Volume With Gaussian CDHR: " << spherical_gaussians_volume(cube3) << "\n";
    std::cout << "Calculated Volume With CRHMC: " << nonspherical_crhmc_volume(cube3) << "\n";
    std::cout << "Expected Volume: " << std::pow(2, 3) << "\n\n";

    HPOLYTOPE cube4 = generate_cube<HPOLYTOPE>(4, false);
    std::cout << "Cube4 \n";
    std::cout << "Calculated Volume With Gaussian CDHR: " << spherical_gaussians_volume(cube4) << "\n";
    std::cout << "Calculated Volume With CRHMC: " << nonspherical_crhmc_volume(cube4) << "\n";
    std::cout << "Expected Volume: " << std::pow(2, 4) << "\n\n";

    HPOLYTOPE skinnycube3 = generate_skinny_cube<HPOLYTOPE>(3, false);
    std::cout << "SkinnyCube3 \n";
    std::cout << "Calculated Volume With Gaussian CDHR: " << spherical_gaussians_volume(skinnycube3) << "\n";
    std::cout << "Calculated Volume With CRHMC: " << nonspherical_crhmc_volume(skinnycube3) << "\n";
    std::cout << "Expected Volume: " << 200 * std::pow(2, 2) << "\n\n";

    HPOLYTOPE P3 = random_hpoly<HPOLYTOPE, PolyRNGType>(3, 12, false);
    std::cout << "Random 3D Hpoly \n";
    std::cout << "Calculated Volume With Gaussian CDHR: " << spherical_gaussians_volume(P3) << "\n";
    std::cout << "Calculated Volume With CRHMC: " << nonspherical_crhmc_volume(P3) << "\n";
    std::cout << "Expected Volume: " << "N/A" << "\n\n";
    
    HPOLYTOPE P4 = random_hpoly<HPOLYTOPE, PolyRNGType>(4, 16, false);
    std::cout << "Random 4D Hpoly \n";
    std::cout << "Calculated Volume With Gaussian CDHR: " << spherical_gaussians_volume(P4) << "\n";
    std::cout << "Calculated Volume With CRHMC: " << nonspherical_crhmc_volume(P4) << "\n";
    std::cout << "Expected Volume: " << "N/A" << "\n\n";
    return 0;
}