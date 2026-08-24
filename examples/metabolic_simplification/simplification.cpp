// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

//
// Usage: ./simplification <model_id.json>
//

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "io/bigg_parser.hpp"
#include "preprocess/metabolic/simplification_clarkson.hpp"
#include "preprocess/metabolic/simplify_and_transform.hpp"
#include "lp_oracles/metabolic_polyoracles.hpp"
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "generators/boost_random_number_generator.hpp"
#include <iostream>
#include <string>

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef MetabolicPolytope<Point> Polytope;
typedef HPolytope<Point> Hpolytope;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DenseMT;
typedef typename Polytope::VT VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNG;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " <model.json>" << std::endl;
        return 1;
    }

    Polytope P = parse_from_json<Point>(argv[1]);

    std::cout << "input polytope" << "\n"
              << " reactions     : " << P.getDimension() << "\n"
              << " metabolites   : " << P.getNumEqualities() << "\n"
              << " finite bounds : " << P.getNumFiniteBounds() << "\n" << std::endl;
    
    ClarksonConfig config;
    config.fix_dimensions = true; 
    
    // Simplifies the polytope, removing redundant bounds and fixing pinned variables, and
    // transforms it to an H-polytope.
    auto [HP, shift, N, simplified] =
        simplify_and_transform<Point, ClarksonSimplifier<Point>, ClarksonConfig>(P, config);

    if (!simplified) {
        std::cerr << "simplification failed" << std::endl;
        return 1;
    }

    std::cout << "\ntransformed polytope\n"
              << " dimension                   : " << HP.dimension() << "\n"
              << " constraints (finite bounds) : " << HP.num_of_hyperplanes() << std::endl;
   
    // Estimates the volume of the original polytope
    unsigned walk_len = 10+HP.dimension()/10;
    NT epsilon = 0.1;

    auto volume = volume_cooling_balls<BallWalk, RNG, Hpolytope>(
        HP, epsilon, walk_len
    ).second;

    std::cout << "\nvolume\n" << "estimate : " << volume << std::endl;
    
    return 0;
}
