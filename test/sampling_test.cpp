// VolEsti (volume computation and sampling library)

// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis
// Copyright (c) 2021 Maios Papachristou

// Licensed under GNU LGPL.3, see LICENCE file


// Edited by HZ on 11.06.2020 - mute doctest.h
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
#include "known_polytope_generators.h"
#include "sampling/sampling.hpp"

#include "diagnostics/univariate_psrf.hpp"


template
<
    typename MT,
    typename WalkType,
    typename Polytope
>
MT get_samples(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    uniform_sampling<WalkType>(randPoints, P, rng, walkL, numpoints,
                               StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}

template
<
    typename MT,
    typename WalkType,
    typename Polytope
>
MT get_samples_boundary(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    uniform_sampling_boundary<WalkType>(randPoints, P, rng, walkL, numpoints,
                                        StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}

template
<
    typename MT,
    typename WalkType,
    typename Polytope
>
MT get_samples_gaussian(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;
    NT a = 1;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    gaussian_sampling<WalkType>(randPoints, P, rng, walkL, numpoints, a,
                               StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}

template <typename NT, typename WalkType = DikinWalk>
void call_test_dikin(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Dikin Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = JohnWalk>
void call_test_john(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing John Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = VaidyaWalk>
void call_test_vaidya(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Vaidya Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = BRDHRWalk>
void call_test_brdhr(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Boundary RDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_boundary<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = BCDHRWalk>
void call_test_bcdhr(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Boundary CDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_boundary<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = GaussianRDHRWalk>
void call_test_grdhr(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Gaussian RDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_gaussian<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = GaussianBallWalk>
void call_test_gbaw(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Gaussian Ball Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_gaussian<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = GaussianHamiltonianMonteCarloExactWalk>
void call_test_ghmc(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Gaussian exact HMC for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_gaussian<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

TEST_CASE("dikin") {
    call_test_dikin<double>();
}

TEST_CASE("john") {
    call_test_john<double>();
}

TEST_CASE("vaidya") {
    call_test_vaidya<double>();
}

TEST_CASE("brdhr") {
    call_test_brdhr<double>();
}

TEST_CASE("bcdhr") {
    call_test_bcdhr<double>();
}

TEST_CASE("grdhr") {
    call_test_grdhr<double>();
}

TEST_CASE("gbaw") {
    call_test_gbaw<double>();
}

TEST_CASE("ghmc") {
    call_test_ghmc<double>();
}
