// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Konstantinos Pallikaris, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SAMPLERS_RANDOM_POINT_GENERATORS_MULTITHREAD_HPP
#define SAMPLERS_RANDOM_POINT_GENERATORS_MULTITHREAD_HPP

#include <iostream>
#include <omp.h>
#include <unistd.h>

template
<
    typename Walk
>
struct RandomPointGeneratorMultiThread
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
                      unsigned int const& num_threads,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        typedef typename Walk::template thread_parameters
        <
            NT,
            Point
        > _thread_parameters;


        omp_set_num_threads(num_threads);
        unsigned int jj = 0, d= P.dimension(), m = P.num_of_hyperplanes;
        std::vector<int> num_points_per_thread(num_threads, 0);

        while (jj < rnum) 
        {
            for (unsigned int i = 0; i < num_threads; i++)
            {
                num_points_per_thread[i]++;
                jj++;
            }
        }

        _thread_parameters thread_random_walk_parameters_temp(d, m);
        Walk walk(P, thread_random_walk_parameters_temp, rng, parameters);

        #pragma omp parallel
        {
            int thread_index = omp_get_thread_num();
            _thread_parameters thread_random_walk_parameters(d, m);
            thread_random_walk_parameters = thread_random_walk_parameters_temp;

            for (unsigned int it = 0; it < num_points_per_thread[thread_index]; it++)
            {
                walk.template apply(P, thread_random_walk_parameters, walk_length, rng);
                policy.apply(randPoints, thread_random_walk_parameters.p);
            }
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
       typedef typename Walk::template thread_parameters
        <
            NT,
            Point
        > _thread_parameters;


        omp_set_num_threads(num_threads);
        unsigned int jj = 0, d= P.dimension(), m = P.num_of_hyperplanes;
        std::vector<int> num_points_per_thread(num_threads, 0);

        while (jj < rnum) 
        {
            for (unsigned int i = 0; i < num_threads; i++)
            {
                num_points_per_thread[i]++;
                jj++;
            }
        }

        _thread_parameters thread_random_walk_parameters_temp(d, m);
        Walk walk(P, thread_random_walk_parameters_temp, rng);

        #pragma omp parallel
        {
            int thread_index = omp_get_thread_num();
            _thread_parameters thread_random_walk_parameters(d, m);
            thread_random_walk_parameters = thread_random_walk_parameters_temp;

            for (unsigned int it = 0; it < num_points_per_thread[thread_index]; it++)
            {
                walk.template apply(P, thread_random_walk_parameters, walk_length, rng);
                policy.apply(randPoints, thread_random_walk_parameters.p);
            }
        }
    }
};


#endif // SAMPLERS_RANDOM_POINT_GENERATORS_MULTITHREAD_HPP
