// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "misc/misc.h"
#include "random_walks/shake_and_bake_walk.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "generators/known_polytope_generators.h"
#include "sampling/sampling.hpp"

#include "diagnostics/univariate_psrf.hpp"

#include "preprocess/feasible_point.hpp"

template
<
    typename MT,
    typename WalkType,
    typename Polytope
>

MT get_samples_shake_and_bake(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);

    auto [boundary_pt, facet_idx] = compute_boundary_point<Point>(P, rng, 1e-7);

    std::list<Point> randPoints;

    shakeandbake_sampling<WalkType>(randPoints,P,rng,walkL,numpoints,boundary_pt, nburns,facet_idx);        

    MT samples(d, numpoints);

    unsigned int jj = 0;
    for (const Point& q : randPoints)
    {
        samples.col(jj++) = q.getCoefficients();
    }
        
    return samples;
}


template <typename NT, typename WalkType = ShakeAndBakeWalk>
void call_test_shake_and_bake(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Running Shake and Bake for H-cube 10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_shake_and_bake<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

TEST_CASE("shake_and_bake") {
    call_test_shake_and_bake<double>();
}
