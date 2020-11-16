// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SAMPLERS_RANDOM_POINT_GENERATORS_HPP
#define SAMPLERS_RANDOM_POINT_GENERATORS_HPP

#include "diagnostics/effective_sample_size.hpp"

template
<
    typename Walk
>
struct RandomPointEfficientGenerator
{
    template
    <
        typename Polytope,
        typename VT,
        typename MT,
        typename RandomNumberGenerator,
        typename Parameters
    >
    static void apply(Polytope& P,
                      VT &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      unsigned int &min_skip,
                      unsigned int &window,
                      MT &randPoints,
                      MT &winPoints,
                      unsigned int const& nburns,
                      bool rounding_requested,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        typedef double NT;
        Walk walk(P, p, rng, parameters);
        int num_points = 100 + 2*int( std::sqrt(NT(P.dimension())) );
        walk.template burnin(P, p, num_points, walk_length, rng);
        bool done = false;
        unsigned int pointer = 0, total_samples = 0;
        int min_eff_samples, min_window = window;
        
        while (!done) 
        {
            for (int i = 0; i < min_window; i++)
            {
                walk.template apply(P, walk_length, rng);
                winPoints.col(pointer + i) = walk.template get_curr_sample();
                if ((i+1)%100 == 0) {
                    std::cout<<"number of sample points = "<<i<<std::endl;
                }
            }
            min_eff_samples = int(ess_univariate_fft<NT, VT>(randPoints).minCoeff());
            std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            if (min_eff_samples <= 0) {
                std::cout<<"[Complete] min_eff_samples = "<<min_eff_samples<<std::endl;
                winPoints.conservativeResize(P.dimension(), winPoints.cols() + window);
                pointer += window;
                //return;
            } else {
                total_samples += min_eff_samples;
                if (winPoints.cols() / min_eff_samples = )
                for (int i = 0; i < ; i++)
                {
                    /* code */
                }
                
            }
            std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            num_samples = NT(total_samples) * (NT(rnum) / NT(min_eff_samples)) + 200 - total_samples;
            std::cout<<"num_samples = "<<num_samples<<std::endl;
            randPoints.conservativeResize(P.dimension(), randPoints.cols() + num_samples);
        }
    }

    template
    <
            typename Polytope,
            typename VT,
            typename MT,
            typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      VT &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      MT &randPoints,
                      unsigned int const& nburns,
                      RandomNumberGenerator &rng)
    {
        typedef double NT;
        Walk walk(P, p, rng);
        bool done = false;
        unsigned int pointer = 0, i = 0, num_samples = rnum, total_samples = 0;
        int min_eff_samples;
        
        while (!done) 
        {
            for (int i = 0; i < num_samples; i++)
            {
                walk.template apply(P, walk_length, rng);
                randPoints.col(pointer + i) = walk.template get_curr_sample();
            }
            total_samples += num_samples;
            pointer = total_samples;
            std::cout<<"total_samples = "<<total_samples<<std::endl;
            min_eff_samples = ess_univariate_fft<NT, VT>(randPoints).minCoeff();
            if (min_eff_samples >= rnum) {
                std::cout<<"[Complete] min_eff_samples = "<<min_eff_samples<<std::endl;
                return;
            }
            std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            num_samples = NT(total_samples) * (NT(rnum) / NT(min_eff_samples)) + 200 - total_samples;
            std::cout<<"num_samples = "<<num_samples<<std::endl;
            randPoints.conservativeResize(P.dimension(), randPoints.cols() + num_samples);
        }
    }
};


template
<
    typename Walk
>
struct RandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator,
        typename Parameters
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, rng, parameters);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, rng);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};

template
<
    typename Walk
>
struct GaussianRandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename NT,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, a_i, rng);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }

    template
    <
            typename Polytope,
            typename Point,
            typename NT,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator,
            typename Parameters
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        Walk walk(P, p, a_i, rng, parameters);

        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};



template <typename Walk>
struct BoundaryRandomPointGenerator
{
    template
    <
            typename Polytope,
            typename Point,
            typename PointList,
            typename WalkPolicy,
            typename RandomNumberGenerator
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, rng);
        Point p1(P.dimension()), p2(P.dimension());
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p1, p2, walk_length, rng);
            policy.apply(randPoints, p1);
            policy.apply(randPoints, p2);
        }
    }
};



#endif // SAMPLERS_RANDOM_POINT_GENERATORS_HPP
