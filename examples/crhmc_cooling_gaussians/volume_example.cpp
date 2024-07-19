// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis
// Copyright (c) 2024 Vladimir Necula

// Contributed and/or modified by Vladimir Necula, as part of Google Summer of
// Code 2024 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "generators/known_polytope_generators.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_gaussians_crhmc.hpp"

#include <iostream>
#include <fstream>
#include "misc/misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt11213b, double> RandomNumberGenerator;
typedef HPolytope<Point> HPOLYTOPE;

void calculateAndVerifyVolume(HPOLYTOPE& polytope) {
    int walk_len = 10;
    NT e = 0.1;

    RandomNumberGenerator rng(polytope.dimension());

    NT volume = volume_cooling_gaussians<HPOLYTOPE, RandomNumberGenerator>(polytope, rng, e, walk_len);

    std::cout << "Volume " << volume << std::endl;
}

int main() {
    
    HPOLYTOPE simplex = generate_simplex<HPOLYTOPE>(2, false);
    std::cout << std::endl << "Simplex: " << std::endl;
    simplex.print();
    calculateAndVerifyVolume(simplex);
    
    HPOLYTOPE cube = generate_cube<HPOLYTOPE>(3, false);
    std::cout << std::endl << "Cube: " << std::endl;
    cube.print();
    calculateAndVerifyVolume(cube);

    HPOLYTOPE cross = generate_cross<HPOLYTOPE>(3, false);
    std::cout << std::endl << "Cross: " << std::endl;
    cross.print();
    calculateAndVerifyVolume(cross);
 
    HPOLYTOPE birkhoff = generate_birkhoff<HPOLYTOPE>(3);
    std::cout << std::endl << "Birkhoff: " << std::endl;
    birkhoff.print();
    calculateAndVerifyVolume(birkhoff);
    
    return 0;
}
