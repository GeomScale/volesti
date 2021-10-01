// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

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
#include "known_polytope_generators.h"

#include "v_polytopes_generators.h"

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
            CHECK(std::abs((volume - expected)/expected) < 0.3);
}

template <class Polytope>
void test_volume(Polytope &P1, Polytope &P2,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& expectedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    // Setup the parameters
    int walk_len = 2;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 105> RNGType;
    typedef IntersectionOfVpoly< Polytope, RNGType > VpIntVp;

    unsigned seed = 105;
    //TODO: low accuracy in high dimensions
    VpIntVp P(P1, P2, seed);
    NT volume = volume_cooling_balls<BallWalk, RNGType>(P, e/2.0, walk_len).second;
    test_values(volume, expectedBall, exact);

    Polytope P11(P1.dimension(), P1.get_mat(), P1.get_vec());
    Polytope P21(P2.dimension(), P2.get_mat(), P2.get_vec());
    VpIntVp P111(P11, P21, seed);
    volume = volume_cooling_balls<CDHRWalk, RNGType>(P111, e/2.0, walk_len).second;
    test_values(volume, expectedCDHR, exact);

    Polytope P12(P1.dimension(), P1.get_mat(), P1.get_vec());
    Polytope P22(P2.dimension(), P2.get_mat(), P2.get_vec());
    VpIntVp P222(P12, P22, seed);
    volume = volume_cooling_balls<RDHRWalk, RNGType>(P222, e/2.0, walk_len).second;
    test_values(volume, expectedRDHR, exact);

    Polytope P13(P1.dimension(), P1.get_mat(), P1.get_vec());
    Polytope P23(P2.dimension(), P2.get_mat(), P2.get_vec());
    VpIntVp P3(P13, P23, seed);
    volume = volume_cooling_balls<BilliardWalk, RNGType>(P3, e/2.0, walk_len).second;
    test_values(volume, expectedBilliard, exact);
}

template <typename NT>
void call_test_vpoly_sphere(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType2;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 105> RNGType;
    typedef VPolytope<Point> Vpolytope;
    typedef IntersectionOfVpoly< Vpolytope, RNGType > VpIntVp;


    std::cout << "--- Testing volume VP1-5-10 , VP2-5-12" << std::endl;
    Vpolytope P1 = random_vpoly<Vpolytope, RNGType2 >(5, 10, 127);
    Vpolytope P2 = random_vpoly<Vpolytope, RNGType2 >(5, 12, 211);

    VpIntVp P(P1, P2, 105);
    if (!P.is_feasible()) {
        std::cout<<"Empty set!"<<std::endl;
        return;
    }

    test_volume(P1, P2, 0.0119774, 0.0136065, 0.013812, 0.0119665, 0.0143023);
}

TEST_CASE("random_vpoly_sphere") {
    call_test_vpoly_sphere<double>();
}
