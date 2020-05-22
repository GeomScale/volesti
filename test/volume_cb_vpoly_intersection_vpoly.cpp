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
            CHECK(std::abs((volume - expected)/expected) < 0.2);
}

template <class VPintersection, class Polytope>
void test_volume(Polytope &P1, Polytope &P2,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& expectedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    VPintersection P;

    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 105> RNGType;

    //TODO: low accuracy in high dimensions
    P.init(P1, P2);
    NT volume = volume_cooling_balls<BallWalk, RNGType>(P, e/2.0, walk_len);
    test_values(volume, expectedBall, exact);

    P1.init(P.dimension(), P1.get_mat(), P1.get_vec());
    P2.init(P.dimension(), P2.get_mat(), P2.get_vec());
    P.init(P1, P2);
    volume = volume_cooling_balls<CDHRWalk, RNGType>(P, e, walk_len);
    test_values(volume, expectedCDHR, exact);

    P1.init(P.dimension(), P1.get_mat(), P1.get_vec());
    P2.init(P.dimension(), P2.get_mat(), P2.get_vec());
    P.init(P1, P2);
    volume = volume_cooling_balls<RDHRWalk, RNGType>(P, e, walk_len);
    test_values(volume, expectedRDHR, exact);

    P1.init(P.dimension(), P1.get_mat(), P1.get_vec());
    P2.init(P.dimension(), P2.get_mat(), P2.get_vec());
    P.init(P1, P2);
    volume = volume_cooling_balls<BilliardWalk, RNGType>(P, e, walk_len);
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
    Vpolytope P1, P2;

    std::cout << "--- Testing volume VP1-5-10 , VP2-5-12" << std::endl;
    P1 = random_vpoly<Vpolytope, RNGType2 >(5, 10, 127);
    P2 = random_vpoly<Vpolytope, RNGType2 >(5, 12, 211);

    VpIntVp P;
    P.init(P1, P2);
    if (!P.is_feasible()) {
        std::cout<<"Empty set!"<<std::endl;
        return;
    }

    test_volume<VpIntVp>(P1, P2, 0.0143023, 0.0143023, 0.0143023, 0.0143023, 0.0143023);
}



TEST_CASE("random_vpoly_sphere") {
    call_test_vpoly_sphere<double>();
}
