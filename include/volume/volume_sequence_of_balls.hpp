// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_SEQUENCE_OF_BALLS_HPP
#define VOLUME_SEQUENCE_OF_BALLS_HPP

#include <iterator>
#include <vector>
#include <list>
#include <math.h>
#include <chrono>

#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/hpolytope.h"
#ifndef DISABLE_LPSOLVE
    #include "convex_bodies/vpolytope.h"
    #include "convex_bodies/zpolytope.h"
    #include "convex_bodies/zonoIntersecthpoly.h"
    #include "convex_bodies/vpolyintersectvpoly.h"
#endif
#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/sampling_policies.hpp"


////////////////////////////// Algorithms


// ----- VOLUME ------ //

template
<
    typename WalkTypePolicy,
    typename Polytope,
    typename RandomNumberGenerator

>
double volume_sequence_of_balls(Polytope& Pin,
                                RandomNumberGenerator &rng,
                                double const& error = 1.0,
                                unsigned int const& walk_length = 1,
                                unsigned int const& n_threads = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    typedef typename Point::FT NT;
    typedef Ball<Point> Ball;
    typedef BallIntersectPolytope<Polytope,Ball> BallPoly;

    typedef typename WalkTypePolicy::template Walk
                                              <
                                                Polytope,
                                                RandomNumberGenerator
                                              > walk;

    typedef RandomPointGenerator<walk> RandomPointGenerator;

    auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    unsigned int rnum = std::pow(error, -2) * 400 * n * std::log(n);
    //RandomNumberGenerator rng(P.dimension());

    //Compute the Chebychev ball (largest inscribed ball) with center and radius
    auto InnerBall = P.ComputeInnerBall();
    if (InnerBall.second < 0.0) return -1.0;

    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(c.getCoefficients());
    c = Point(n);

    // Scale by number of threads and prevent edge case rnum=0 from producing overflow later
    rnum = rnum >= 2*n_threads ? rnum/n_threads : 2u;
    NT vol = NT(0);

    // Perform the procedure for a number of threads and then take the average
    for (auto t=0u; t<n_threads; t++)
    {
        // Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
#ifdef VOLESTI_DEBUG
        std::cout<<"\nGenerate the first random point in P"<<std::endl;
#endif
        Point p = GetPointInDsphere<Point>::apply(P.dimension(), radius, rng);
        std::list<Point> randPoints; //ds for storing rand points

        PushBackWalkPolicy push_back_policy;
        RandomPointGenerator::apply(P, p, 1, 50*n, randPoints, push_back_policy, rng);

#ifdef VOLESTI_DEBUG
        double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
#endif
        RandomPointGenerator::apply(P, p, rnum-1, walk_length, randPoints,
                                    push_back_policy, rng);

#ifdef VOLESTI_DEBUG
        double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "First random points construction time = "
                  << tstop2 - tstart2 << std::endl;
#endif

        // Construct the sequence of balls
        // a. compute the radius of the largest ball
        NT current_dist, max_dist=NT(0);
        for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit)
        {
            current_dist = (*pit).squared_length();
            if (current_dist > max_dist)
            {
                max_dist=current_dist;
            }
        }
        max_dist = std::sqrt(max_dist);
#ifdef VOLESTI_DEBUG
        std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist
                <<std::endl;
        std::cout<<"\nConstructing the sequence of balls"<<std::endl;
        std::cout<<"---------"<<std::endl;
#endif

        //
        // b. Number of balls
        int nb1 = n * (std::log(radius)/std::log(2.0));
        int nb2 = std::ceil(n * (std::log(max_dist)/std::log(2.0)));

        std::vector<Ball> balls;

        for (auto i=nb1; i<=nb2; ++i)
        {
            if (i == nb1)
            {
                balls.push_back(Ball(c,radius*radius));
                vol = (std::pow(M_PI,n/2.0)*(std::pow(balls[0].radius(), n) ) )
                       / (tgamma(n/2.0+1));
            } else {
                balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
            }
        }
        assert(!balls.empty());

        // Estimate Vol(P)
        typename std::vector<Ball>::iterator bit2=balls.end();
        bit2--;

        while (bit2!=balls.begin())
        {
            //each step starts with some random points in PBLarge stored
            //in list "randPoints", these points have been generated in a
            //previous step

            BallPoly PBLarge(P,*bit2);
            --bit2;
            BallPoly PBSmall(P,*bit2);

#ifdef VOLESTI_DEBUG
            std::cout<<"("<<balls.end()-bit2<<"/"<<balls.end()-balls.begin()
                     <<")Ball ratio radius="
                     <<PBLarge.second().radius()<<","
                     <<PBSmall.second().radius()<<std::endl;
            std::cout<<"Points in PBLarge="<<randPoints.size()<<std::endl;
#endif

            // choose a point in PBLarge to be used to generate more rand points
            Point p_gen = *randPoints.begin();

            // num of points in PBSmall and PBLarge
            unsigned int nump_PBSmall = 0;
            unsigned int nump_PBLarge = randPoints.size();

            //keep the points in randPoints that fall in PBSmall
            typename std::list<Point>::iterator rpit=randPoints.begin();
            while (rpit!=randPoints.end())
            {
                if (PBSmall.second().is_in(*rpit) == 0)//not in
                {
                    rpit=randPoints.erase(rpit);
                } else {
                    ++nump_PBSmall;
                    ++rpit;
                }
            }

#ifdef VOLESTI_DEBUG
            std::cout<<"Points in PBSmall="<<randPoints.size()
                     <<"\nRatio= "<<NT(nump_PBLarge)/NT(nump_PBSmall)
                     <<std::endl;
            std::cout<<"Generate "<<rnum-nump_PBLarge<<" more "<<std::endl;
#endif

            CountingWalkPolicy<BallPoly> counting_policy(nump_PBSmall, PBSmall);
            RandomPointGenerator::apply(PBLarge, p_gen, rnum-nump_PBLarge,
                                        walk_length, randPoints,
                                        counting_policy, rng);

            nump_PBSmall = counting_policy.get_nump_PBSmall();

            vol *= NT(rnum)/NT(nump_PBSmall);

#ifdef VOLESTI_DEBUG
            std::cout<<nump_PBSmall<<"/"<<rnum<<" = "<<NT(rnum)/nump_PBSmall
                     <<"\ncurrent_vol = "<<vol
                     <<"\n--------------------------"<<std::endl;
#endif
            //don't continue in pairs of balls that are almost inside P, i.e. ratio ~= 2
        }
    }
#ifdef VOLESTI_DEBUG
    std::cout<<"rand points = "<<rnum<<std::endl;
    std::cout<<"walk len = "<<walk_length<<std::endl;
    std::cout<<"volume computed: "<<vol<<std::endl;
#endif

    return vol;
}


template
<
    typename WalkTypePolicy = CDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b, double>,
    typename Polytope
>
double volume_sequence_of_balls(Polytope &Pin,
                                double const& error = 1.0,
                                unsigned int const& walk_length = 1,
                                unsigned int const& n_threads = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_sequence_of_balls<WalkTypePolicy>(Pin, rng, error,
                                                    walk_length, n_threads);
}


template
<
    typename WalkTypePolicy = CDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b, double>,
    typename Polytope
>
double volume_sequence_of_balls(Polytope &Pin,
                                Cartesian<double>::Point const& interior_point,
                                unsigned int const& walk_length = 1,
                                double const& error = 1.0,
                                unsigned int const& n_threads = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    Pin.set_interior_point(interior_point);

    return volume_sequence_of_balls<WalkTypePolicy>(Pin, rng, error,
                                                    walk_length, n_threads);
}

#endif // VOLUME_SEQUENCE_OF_BALLS_HPP
