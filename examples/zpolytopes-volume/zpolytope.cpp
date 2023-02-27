// VolEsti ( volume computation and sampling library)

// Copyright (c) 2012-2023 Vissarion Fisikopoulos
// Copyright (c) 2023 Apostolos Chalkis

// Conributed and/or modified by Sarthak Bhandari
// Licensed under GNU LGPL.3, see LICENCE file


#include "Eigen/Eigen"
#include <vector>
#include "cartesian_geom/cartesian_kernel.h"
#include "zpolytope/zpolytope.h"
#include "known_zpolytope_generators.h"
#include "random_walks/random_walks.hpp"
#include "misc.h"

#include "volume_cooling_gaussians.hpp"
#include "volume_cooling_balls.hpp"

typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

void calculateVolumes(const ZPolytope<Kernel> &ZP) {
    // Setup parameters for calculating volume
    int walk_len = 10 + ZP.dimension()/10;
    NT e = 0.1;

    // Calculating volume of the passed polytope
    NT volume1 = volume_cooling_balls<BallWalk, RNGType, ZPolytope<Kernel>>(ZP, e, walk_len).second;
    NT volume2 = volume_cooling_gaussians<GaussianBallWalk, RNGType, ZPolytope<Kernel>>(ZP, e, walk_len);

    std::cout << "\t Using Cooling Balls method: " << volume1 << "\n";
    std::cout << "\t Using Cooling Gaussians method: " << volume2 << "\n";
}

int main(int argc, char* argv[]) {
    // Generating a 4-dimensional ZPolytope centered at origin
    ZPolytope<Kernel> ZP1 = generate_hypercube<ZPolytope<Kernel>>(4);
    std::cout << "ZPolytope ZP1: \n";
    ZP1.print();
    std::cout << "\n";

    std::cout << "Volume of ZP1: \n";
    calculateVolumes(ZP1);

    // Generating a 3-dimensional ZPolytope with random generator matrix
    int dim = 3;
    Eigen::MatrixXd M(dim, dim);
    M.setRandom();
    ZPolytope<Kernel> ZP2(M);
    std::cout << "ZPolytope ZP2: \n";
    ZP2.print();
    std::cout << "\n";

    std::cout << "Volume of ZP2: \n";
    calculateVolumes(ZP2);

    return 0;
}
