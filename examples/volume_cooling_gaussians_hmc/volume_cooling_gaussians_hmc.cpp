#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"

#include <iostream>
#include <fstream>
#include "misc/misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef HPolytope<Point> HPOLYTOPE;

void verify_volume(HPOLYTOPE& polytope, const std::string& description, NT expectedVolume, NT epsilon, NT L) {
    int walk_len = 1;
    NT e = 0.1;

    std::cout << "Calculating volume for " << description << ":\n";
    NT volume = volume_cooling_gaussians<GaussianHamiltonianMonteCarloExactWalk, RNGType>(polytope, e, walk_len, L);

    std::cout << "Estimated volume: " << volume << "\n";
    std::cout << "Expected volume: " << expectedVolume << "\n";
    NT error = std::abs(volume - expectedVolume);
    
    if (error <= epsilon)
        std::cout << "Result is within the acceptable error margin (" << epsilon << ").\n\n";
    else
        std::cout << "Error (" << error << ") is larger than acceptable margin (" << epsilon << ").\n\n";
}

int main() {
    const NT epsilon = 0.1;
    
    // 3-dimensional cube
    HPOLYTOPE cube = generate_cube<HPOLYTOPE>(3, false);
    verify_volume(cube, "3-dimensional cube", 1.0, epsilon, 4);

    // 3-dimensional cross polytope
    HPOLYTOPE crossPolytope = generate_cross<HPOLYTOPE>(3, false);
    verify_volume(crossPolytope, "3-dimensional cross polytope", 4.0 / 6.0, epsilon, 4);

    // 3-dimensional simplex
    HPOLYTOPE simplex = generate_simplex<HPOLYTOPE>(3, false);
    verify_volume(simplex, "3-dimensional simplex", 1.0 / 6.0, epsilon, 10);

    // 4-dimensional birkhoff  
    HPOLYTOPE birkhoffPolytope = generate_birkhoff<HPOLYTOPE>(3);
    verify_volume(birkhoffPolytope, "Birkhoff polytope (dim 4)", 1.0 / 16.0, epsilon, 2);  // Theoretical volume not easily available

    return 0;
}

