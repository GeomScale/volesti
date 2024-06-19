// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis

// Contributed by Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>

#include <boost/random.hpp>

#include "misc/misc.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"

#include "preprocess/max_inscribed_ball.hpp"

#include "generators/known_polytope_generators.h"
#include "generators/h_polytopes_generator.h"

template <typename NT>
void call_test_max_ball() {
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope <Point> Hpolytope;
    typedef boost::mt19937 PolyRNGType;
    Hpolytope P;

    std::cout << "\n--- Testing rounding of H-skinny_cube5" << std::endl;
    P = skinny_random_hpoly<Hpolytope, NT, PolyRNGType>(100, 1000, true, 2000.0);
    
    std::pair<Point, NT> InnerBall = P.ComputeInnerBall();

    auto [center, radius, converged] =  max_inscribed_ball(P.get_mat(), P.get_vec(), 500, 1e-06);
    
    std::cout<<"[1] center: "<<InnerBall.first.getCoefficients().transpose()<<std::endl;
    std::cout<<"[1] radius: "<<InnerBall.second<<std::endl;

    std::cout<<"[2] center: "<<center.transpose()<<std::endl;
    std::cout<<"[2] radius: "<<radius<<std::endl;
    CHECK(false);
}

TEST_CASE("test_max_ball") {
    call_test_max_ball<double>();
}