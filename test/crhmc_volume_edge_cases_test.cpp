// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

// Edge cases of the CRHMC cooling-Gaussians volume algorithm.

#include "doctest.h"
#include <stdexcept>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_cooling_gaussians_crhmc.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef HPolytope<Point> Hpolytope;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

void call_test_crhmc_inner_ball_failure()
{
    // An empty polytope (x <= 1 and x >= 2) has no inner ball: throw
    // instead of returning -1
    MT A(4, 2);
    A << 1, 0, 0, 1, -1, 0, 0, -1;
    VT b(4);
    b << 1, 1, -2, 1;
    Hpolytope P(2, A, b);
    RNGType rng(2);
    CHECK_THROWS_AS((volume_cooling_gaussians<Hpolytope, RNGType>(P, rng, 0.1, 1)),
                    std::runtime_error);
}

TEST_CASE("crhmc_inner_ball_failure") {
    call_test_crhmc_inner_ball_failure();
}
