// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include "new_volume.hpp"
#include "new_gaussian_volume.hpp"
#include "new_cooling_balls.hpp"
#include "new_cooling_hpoly.hpp"
#include "exact_vols.h"
#include "z_polytopes_gen.h"

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
            CHECK(std::abs((volume - exact)/exact) < 0.12);
}

template <class Polytope>
void test_volume(Polytope &P,
                 double const& expectedBall,
                 double const& expectedCDHR,
                 double const& expectedRDHR,
                 double const& expectedBilliard,
                 double const& exact)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef HPolytope<Point> Hpolytope;

    // Setup the parameters
    int walk_len = 1;
    NT e = 0.1, volume;

    // Estimate the volume
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 5> RNGType;

    //TODO: low accuracy in high dimensions
    //NT volume = volume_cooling_balls<BallWalk, RNGType>(HP, e, walk_len);
    //test_values(volume, expectedBall, exact);
    P.init(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_hpoly<CDHRWalk, RNGType, Hpolytope>(P, e, walk_len);
    test_values(volume, expectedCDHR, exact);

    P.init(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_hpoly<RDHRWalk, RNGType, Hpolytope>(P, e, walk_len);
    test_values(volume, expectedRDHR, exact);

    P.init(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_hpoly<BilliardWalk, RNGType, Hpolytope>(P, e, walk_len);
    test_values(volume, expectedBilliard, exact);

    //cooling balls//
    P.init(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_balls<CDHRWalk, RNGType>(P, e, walk_len);
    test_values(volume, expectedCDHR, exact);

    P.init(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_balls<RDHRWalk, RNGType>(P, e, walk_len);
    test_values(volume, expectedRDHR, exact);

    P.init(P.dimension(), P.get_mat(), P.get_vec());
    volume = volume_cooling_balls<BilliardWalk, RNGType>(P, e, walk_len);
    test_values(volume, expectedBilliard, exact);
}

template <typename NT>
void call_test_uniform_generator(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937            RNGType;
    typedef Zonotope<Point> zonotope;
    zonotope P;

    P = gen_zonotope_uniform<zonotope, RNGType>(5, 10, 127);
    NT exact_vol = exact_zonotope_vol<NT>(P);
    test_volume(P, exact_vol, exact_vol, exact_vol, exact_vol, exact_vol);

    P = gen_zonotope_uniform<zonotope, RNGType>(10, 15, 211);
    exact_vol = exact_zonotope_vol<NT>(P);
    test_volume(P, exact_vol, exact_vol, exact_vol, exact_vol, exact_vol);
}


TEST_CASE("uniform_zonotopes") {
    call_test_uniform_generator<double>();
    //call_test_cube_float<float>();
}

