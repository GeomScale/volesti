// VolEsti (volume computation and sampling library)

// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis
// Copyright (c) 2022 Konstantinos Pallikaris

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
#include "known_polytope_generators.h"
#include "sampling/random_point_generators_multithread.hpp"

#include "diagnostics/univariate_psrf.hpp"


template <typename MT,
          typename WalkTypePolicy,
          typename Point,
          typename Polytope,
          typename RandomNumberGenerator
         >
MT get_uniform_samples(Polytope &P,
                         RandomNumberGenerator &rng,
                         const unsigned int &walk_len,
                         const unsigned int &N)
{
    typedef typename Point::FT NT;

    typedef typename WalkTypePolicy::template thread_parameters
        <
            NT,
            Point
        > _thread_parameters;
    
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator,
                    _thread_parameters
            > walk;

    PushBackWalkPolicy push_back_policy;

    Point p = P.ComputeInnerBall().first;
    std::list<Point> randPoints;
    unsigned int num_threads = 2;

    typedef RandomPointGeneratorMultiThread <walk> _RandomPointGeneratorMultiThread;
    
    _RandomPointGeneratorMultiThread::apply(P, p, N, walk_len, num_threads, randPoints,
                                            push_back_policy, rng);

    unsigned int d = P.dimension();
    
    MT samples(d, N);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}


template <typename MT,
          typename WalkTypePolicy,
          typename Point,
          typename Polytope,
          typename NT,
          typename RandomNumberGenerator
         >
MT get_gaussian_samples(Polytope &P,
                          NT const& a_i,
                          RandomNumberGenerator &rng,
                          const unsigned int &walk_len,
                          const unsigned int &N)
{
    typedef typename WalkTypePolicy::template thread_parameters
        <
            NT,
            Point
        > _thread_parameters;
    
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator,
                    _thread_parameters
            > walk;

    PushBackWalkPolicy push_back_policy;

    Point p = P.ComputeInnerBall().first;
    std::list<Point> randPoints;
    unsigned int num_threads = 2;

    typedef GaussianPointGeneratorMultiThread <walk> _GaussianPointGeneratorMultiThread;
    
    _GaussianPointGeneratorMultiThread::apply(P, p, a_i, N, walk_len, num_threads, randPoints,
                                              push_back_policy, rng);

    unsigned int d = P.dimension();

    MT samples(d, N);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}


template <typename NT, typename WalkType>
void call_test_parallel_uniform_random_walk(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    Hpolytope P;
    unsigned int d = 10, walk_len = 5, N = 5000;

    std::cout << "--- Testing Boundary RDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    //P.ComputeInnerBall();
    RNGType rng(P.dimension());

    MT samples = get_uniform_samples<MT, WalkType, Point>(P, rng, walk_len, N);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType>
void call_test_parallel_gaussian_random_walk(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    Hpolytope P;
    unsigned int d = 10, walk_len = 5, N = 5000;
    NT a_i = 1;

    std::cout << "--- Testing Boundary RDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    //P.ComputeInnerBall();
    RNGType rng(P.dimension());

    MT samples = get_gaussian_samples<MT, WalkType, Point>(P, a_i, rng, walk_len, N);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}


TEST_CASE("multithread_brdhr") {
    call_test_parallel_uniform_random_walk<double, BRDHRWalk_multithread>();
}

TEST_CASE("multithread_bcdhr") {
    call_test_parallel_uniform_random_walk<double, BCDHRWalk_multithread>();
}

TEST_CASE("multithread_gcdhr") {
    call_test_parallel_gaussian_random_walk<double, GaussianCDHRWalk_multithread>();
}

TEST_CASE("multithread_grdhr") {
    call_test_parallel_gaussian_random_walk<double, GaussianRDHRWalk_multithread>();
}

TEST_CASE("multithread_biw") {
    call_test_parallel_uniform_random_walk<double, BilliardWalk_multithread>();
}

TEST_CASE("multithread_cdhr") {
    call_test_parallel_uniform_random_walk<double, CDHRWalk_multithread>();
}

TEST_CASE("multithread_rdhr") {
    call_test_parallel_uniform_random_walk<double, RDHRWalk_multithread>();
}
