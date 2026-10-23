// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

// Headers used on their own, without random_walks/random_walks.hpp, as the
// bindings do to keep their include set small. Built with DISABLE_LPSOLVE,
// as the Octave binding is.

#include "doctest.h"
#include <cmath>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
#include "random_walks/uniform_accelerated_billiard_walk.hpp"
#include "volume/volume_cooling_balls.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef HPolytope<Point> Hpolytope;

void call_test_accelerated_billiard()
{
    // uniform_accelerated_billiard_walk.hpp used compute_diameter, and
    // volume_cooling_balls.hpp PushBackWalkPolicy, without including them
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    MT A(6, 3);
    A << MT::Identity(3, 3), -MT::Identity(3, 3);
    VT b = VT::Ones(6);
    Hpolytope P(3, A, b);
    NT volume = volume_cooling_balls<AcceleratedBilliardWalk, RNGType>(P, 0.1, 1).second;
    std::cout << "Computed volume " << volume << ", exact 8" << std::endl;
    CHECK(std::abs(volume - 8) / 8 < 0.2);
}

TEST_CASE("accelerated_billiard") {
    call_test_accelerated_billiard();
}
