// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_COOLING_GAUSSIANS_GENERAL_MASS_MATRIX_HPP
#define VOLUME_COOLING_GAUSSIANS_GENERAL_MASS_MATRIX_HPP

//#define VOLESTI_DEBUG

#include <iterator>
#include <vector>
#include <list>
#include <math.h>
#include <chrono>

#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/gaussian_helpers.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "random_walks/gaussian_hamiltonian_monte_carlo_exact_walk.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/math_helpers.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "Eigen/Eigen"



////////////////////////////// Algorithms

// Gaussian Anealling

// Implementation is based on algorithm from paper "A practical volume algorithm",
// Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society 2015
// Ben Cousins, Santosh Vempala

// Compute the sequence of spherical gaussians
template
<
    typename WalkType,
    typename RandomPointGenerator,
    typename Polytope,
    typename NT,
    typename RandomNumberGenerator
>
void compute_general_annealing_schedule(Polytope const& P,
                                NT const& ratio,
                                NT const& C,
                                NT const& frac,
                                unsigned int const& N,
                                unsigned int const& walk_length,
                                NT const& chebychev_radius,
                                NT const& error,
                                std::vector<NT>& a_vals,
                                RandomNumberGenerator& rng,
				Eigen::MatrixXd& covar)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::VT VT;

    // Compute the first gaussian
    get_first_gaussian(P, frac, chebychev_radius, error, a_vals);

#ifdef VOLESTI_DEBUG
    std::cout<<"first gaussian computed\n"<<std::endl;
#endif

    NT a_stop = 0.0;
    const NT tol = 0.001;
    unsigned int it = 0;
    unsigned int n = P.dimension();
    const unsigned int totalSteps = ((int)150/((1.0 - frac) * error))+1;

    if (a_vals[0]<a_stop) a_vals[0] = a_stop;

#ifdef VOLESTI_DEBUG
    std::cout<<"Computing the sequence of gaussians..\n"<<std::endl;
#endif

    Point p(n);
    std::vector<Point> points;

    while (true)
    {
        // Compute the next gaussian
        NT next_a = get_next_gaussian<RandomPointGenerator>
                      (P, p, a_vals[it], N, ratio, C, walk_length, rng);

        NT curr_fn = 0;
        NT curr_its = 0;
        auto steps = totalSteps;

        WalkType walk(P, p, a_vals[it], rng);
        //TODO: test update delta here?

        update_delta<WalkType>
                ::apply(walk, 4.0 * chebychev_radius
                        / std::sqrt(std::max(NT(1.0), a_vals[it]) * NT(n)));
	points.push_back(p);

        // Compute some ratios to decide if this is the last gaussian
        for (unsigned  int j = 0; j < steps; j++)
        {
            walk.template apply(P, p, a_vals[it], walk_length, rng);
            curr_its += 1.0;
            curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
        }

        // Remove the last gaussian.
        // Set the last a_i equal to zero
        if (next_a>0 && curr_fn/curr_its>(1.0+tol))
        {
            a_vals.push_back(next_a);
            it++;
        } else if (next_a <= 0)
        {
            a_vals.push_back(a_stop);
            it++;
            break;
        } else {
            a_vals[it] = a_stop;
            break;
        }
    }
    Point mean_point(n);
    for(int i=0;i<n;i++)
    {
	    mean_point.set_coord(i, 0);
	    NT sum = 0;
	    for(int j=0;j<points.size();j++)
		    sum += points[j][i];
	    mean_point.set_coord(i, sum/n);
    }
    for(int i=0;i<n;i++)
    {
	    for(int j=0;j<n;j++)
	    {
		    covar(i, j) = 0;
		    for(int k=0;k<points.size();k++)
			    covar(i, j) += (points[k][i]-mean_point[i])*(points[k][j]-mean_point[j]);
		    covar(i, j) = covar(i, j)/(n-1);
	    }
    }
}

template
<
    typename WalkTypePolicy,
    typename Polytope,
    typename RandomNumberGenerator

>
double volume_general_cooling_gaussians(Polytope const& Pin,
                                RandomNumberGenerator& rng,
                                double const& error = 0.1,
                                unsigned int const& walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef typename Polytope::VT 	VT;
    typedef typename WalkTypePolicy::template Walk
                                              <
                                                    Polytope,
                                                    RandomNumberGenerator
                                              > WalkType;
    typedef GaussianRandomPointGenerator<WalkType> RandomPointGenerator;

    //const NT maxNT = std::numeric_limits<NT>::max();//1.79769e+308;
    //const NT minNT = std::numeric_limits<NT>::min();//-1.79769e+308;

    auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    unsigned int m = P.num_of_hyperplanes();
    gaussian_annealing_parameters<NT> parameters(P.dimension());
    //RandomNumberGenerator rng(n);

    // Consider Chebychev center as an internal point
    auto InnerBall = P.ComputeInnerBall();
    if (InnerBall.second < 0.0) return -1.0;

    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(c.getCoefficients());

    // Computing the sequence of gaussians
#ifdef VOLESTI_DEBUG
    std::cout<<"\n\nComputing annealing...\n"<<std::endl;
    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
#endif

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = parameters.ratio;
    NT C = parameters.C;
    unsigned int N = parameters.N;
    Eigen::MatrixXd covar(n, n);

    compute_general_annealing_schedule
    <
        WalkType,
        RandomPointGenerator
    >(P, ratio, C, parameters.frac, N, walk_length, radius, error, a_vals, rng, covar);

#ifdef VOLESTI_DEBUG
    std::cout<<"All the variances of schedule_annealing computed in = "
            << (double)clock()/(double)CLOCKS_PER_SEC-tstart2<<" sec"<<std::endl;
    auto j=0;
    for (auto avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++){
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;
#endif

    // Initialization for the approximation of the ratios
    unsigned int W = parameters.W;
    unsigned int mm = a_vals.size()-1;
    std::vector<NT> last_W2(W,0);
    std::vector<NT> fn(mm,0);
    std::vector<NT> its(mm,0);
    VT lamdas;
    lamdas.setZero(m);
    NT vol = std::pow(M_PI/a_vals[0], (NT(n))/2.0);
    Point p(n); // The origin is the Chebychev center of the Polytope
    for(int i=0;i<n;i++)p.set_coord(i, covar(i, i));
    unsigned int i=0;

    typedef typename std::vector<NT>::iterator viterator;
    viterator itsIt = its.begin();
    viterator avalsIt = a_vals.begin();
    viterator minmaxIt;

#ifdef VOLESTI_DEBUG
    std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    std::cout<<"computing ratios..\n"<<std::endl;
#endif

    //iterate over the number of ratios
    for (viterator fnIt = fn.begin();
         fnIt != fn.end();
         fnIt++, itsIt++, avalsIt++, i++)
    {
        //initialize convergence test
        bool done = false;
        NT curr_eps = error/std::sqrt((NT(mm)));
        NT min_val = std::numeric_limits<NT>::min();
        NT max_val = std::numeric_limits<NT>::max();
        unsigned int min_index = W-1;
        unsigned int max_index = W-1;
        unsigned int index = 0;
        unsigned int min_steps = 0;
        std::vector<NT> last_W = last_W2;

        // Set the radius for the ball walk
        WalkType walk(P, p, *avalsIt, rng);

        update_delta<WalkType>
                ::apply(walk, 4.0 * radius
                         / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n)));

        while (!done || (*itsIt)<min_steps)
        {
            walk.template apply(P, p, *avalsIt, walk_length, rng);

            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
            NT val = (*fnIt) / (*itsIt);

            last_W[index] = val;
            if (val <= min_val)
            {
                min_val = val;
                min_index = index;
            } else if (min_index == index)
            {
                minmaxIt = std::min_element(last_W.begin(), last_W.end());
                min_val = *minmaxIt;
                min_index = std::distance(last_W.begin(), minmaxIt);
            }

            if (val >= max_val)
            {
                max_val = val;
                max_index = index;
            } else if (max_index == index)
            {
                minmaxIt = std::max_element(last_W.begin(), last_W.end());
                max_val = *minmaxIt;
                max_index = std::distance(last_W.begin(), minmaxIt);
            }

            if ( (max_val-min_val)/max_val <= curr_eps/2.0 )
            {
                done=true;
            }

            index = index%W + 1;
            if (index == W) index = 0;
        }
#ifdef VOLESTI_DEBUG
        std::cout << "ratio " << i << " = " << (*fnIt) / (*itsIt)
                  << " N_" << i << " = " << *itsIt << std::endl;
#endif
        vol *= ((*fnIt) / (*itsIt));
    }

#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        auto steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
#endif

    return vol;
}


template
<
    typename WalkTypePolicy = GaussianCDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b, double>,
    typename Polytope
>
double volume_general_cooling_gaussians(Polytope &Pin,
                                 double const& error = 0.1,
                                 unsigned int const& walk_length = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_general_cooling_gaussians<WalkTypePolicy>(Pin, rng, error, walk_length);
}


template
<
    typename WalkTypePolicy = GaussianCDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b, double>,
    typename Polytope
>
double volume_general_cooling_gaussians(Polytope &Pin,
                                Cartesian<double>::Point const& interior_point,
                                unsigned int const& walk_length = 1,
                                double const& error = 0.1)
{
    RandomNumberGenerator rng(Pin.dimension());
    Pin.set_interior_point(interior_point);
    
    return volume_general_cooling_gaussians<WalkTypePolicy>(Pin, rng, error, walk_length);
}


#endif // VOLUME_COOLING_GAUSSIANS_GENERAL_MASS_MATRIX_HPP
