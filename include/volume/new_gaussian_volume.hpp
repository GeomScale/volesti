// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NEW_GAUSSIAN_VOLUME_HPP
#define NEW_GAUSSIAN_VOLUME_HPP

#include <iterator>
#include <vector>
#include <list>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "ball.h"
#include "ballintersectconvex.h"
#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
#include "gaussian_samplers.h"
#include "gaussian_annealing.h"

#include "new_volume.hpp"


// Pick a point from the distribution exp(-a_i||x||^2) on the chord
template
<
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
void chord_random_point_generator_exp(Point &lower,
                                      Point & upper,
                                      const NT &a_i,
                                      Point &p,
                                      RandomNumberGenerator& rng)
{
    NT r, r_val, fn;
    const NT tol = 0.00000001;
    Point bef = upper - lower;
    // pick from 1-dimensional gaussian if enough weight is inside polytope P
    if (a_i > tol && std::sqrt(bef.squared_length()) >= (2.0 / std::sqrt(2.0 * a_i)))
    {
        Point a = -1.0 * lower;
        Point b = (1.0 / std::sqrt(bef.squared_length())) * bef;
        Point z = (a.dot(b) * b) + lower;
        NT low_bd = (lower[0] - z[0]) / b[0];
        NT up_bd = (upper[0] - z[0]) / b[0];
        while (true) {
            r = rng.sample_ndist();//rdist(rng2);
            r = r / std::sqrt(2.0 * a_i);
            if (r >= low_bd && r <= up_bd) {
                break;
            }
        }
        p = (r * b) + z;

    // select using rejection sampling from a bounding rectangle
    } else {
        NT M = get_max(lower, upper, a_i);
        while (true) {
            r = rng.sample_urdist();//urdist(rng2);
            Point pef = r * upper;
            p = ((1.0 - r) * lower) + pef;
            r_val = M * rng.sample_urdist();//urdist(var.rng);
            fn = eval_exp(p, a_i);
            if (r_val < fn) {
                break;
            }
        }
    }
}

/////////////////// Random Walks

// ball walk with gaussian target distribution
struct GaussianBallWalk
{

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    Walk (Polytope const&, Point&, RandomNumberGenerator&)
    {}

    Walk (BallPolytope const&, Point &, RandomNumberGenerator &)
    {}

    template<typename BallPolytope>
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        const NT delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());

        for (auto j=0u; j<walk_length; ++j)
        {
            Point y = GetPointInDsphere<Point>::apply(P.dimension(),
                                                      delta,
                                                      rng);
            y += p;
            if (P.is_in(y) == -1)
            {
                NT f_x = eval_exp(p, a_i);
                NT f_y = eval_exp(y, a_i);
                NT rnd = rng.sample_urdist();
                if (rnd <= f_y / f_x) {
                    p = y;
                }
            }
        }
    }
};

};

// random directions hit-and-run walk with uniform target distribution
struct GaussianRDHRWalk
{

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;

    Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng)
    {}

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      NT const& a_i,
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            Point v = GetDirection<Point>::apply(p.dimension(), rng);
            std::pair <NT, NT> dbpair = P.line_intersect(p, v);

            NT min_plus = dbpair.first;
            NT max_minus = dbpair.second;
            Point upper = (min_plus * v) + p;
            Point lower = (max_minus * v) + p;

            chord_random_point_generator_exp(lower, upper, a_i, p, rng);
        }
    }
};

};

// random directions hit-and-run walk with uniform target distribution
struct GaussianCDHRWalk
{

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            auto rand_coord_prev = _rand_coord;
            _rand_coord = rng.sample_uidist();
            NT kapa = rng.sample_urdist();
            std::pair<NT, NT> bpair = P.line_intersect_coord(_p,
                                                             _p_prev,
                                                             _rand_coord,
                                                             rand_coord_prev,
                                                             _lamdas);
            _p_prev = _p;
            _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                         * (bpair.second - bpair.first));
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           Point &p,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _rand_coord = rng.sample_uidist();
        NT kapa = rng.sample_urdist();
        _p = p;
        std::pair<NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord, _lamdas);
        _p_prev = _p;
        _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                    * (bpair.second - bpair.first));
    }

    unsigned int _rand_coord;
    Point _p;
    Point _p_prev;
    typename Point::Coeff _lamdas;
};

};


////////////////////////////// Random Point Generators
///

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
        Walk walk(P, p, rng);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, a_i, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};


////////////////////////////// Algorithms



///////////////////////////////////////////
// Gaussian Anealling


// Implementation is based on algorithm from paper "A practical volume algorithm",
// Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society 2015
// Ben Cousins, Santosh Vempala

//An implementation of Welford's algorithm for mean and variance.
template <typename NT>
std::pair<NT, NT> getMeanVariance2(std::vector<NT>& vec)
{
    NT mean = 0;
    NT M2 = 0;
    NT variance = 0;
    NT delta;

    unsigned int i=0;
    for (auto vecit = vec.begin(); vecit!=vec.end(); vecit++, i++){
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }

    return std::pair<NT, NT> (mean, variance);
}


// Compute the first variance a_0 for the starting gaussian
template <typename Polytope, typename Parameters, typename NT>
void get_first_gaussian(Polytope & P,
                        NT const& frac,
                        Parameters const& var,
                        NT & error,
                        std::vector<NT> & a_vals)
{

    unsigned int i;
    const unsigned int maxiter = 10000;
    typedef typename std::vector<NT>::iterator viterator;
    NT tol;
    // if tol is smaller than 1e-6 no convergence can be obtained when float is used
    if (eqTypes<float, NT>()) {
        tol = 0.001;
    } else {
        tol = 0.0000001;
    }

    NT sum, lower = 0.0, upper = 1.0, mid;
    std::vector <NT> dists = P.get_dists(var.che_rad);

    // Compute an upper bound for a_0
    for (i= 1; i <= maxiter; ++i) {
        sum = 0.0;
        for (viterator it = dists.begin(); it != dists.end(); ++it)
        {
            sum += std::exp(-upper * std::pow(*it, 2.0))
                    / (2.0 * (*it) * std::sqrt(M_PI * upper));
        }

        if (sum > frac * error)
        {
            upper = upper * 10;
        } else {
            break;
        }
    }

    if (i == maxiter) {
#ifdef VOLESTI_DEBUG
        std::cout << "Cannot obtain sharp enough starting Gaussian" << std::endl;
#endif
        return;
    }

    //get a_0 with binary search
    while (upper - lower > tol) {
        mid = (upper + lower) / 2.0;
        sum = 0.0;
        for (viterator it = dists.begin(); it != dists.end(); ++it) {
            sum += std::exp(-mid * std::pow(*it, 2.0))
                    / (2.0 * (*it) * std::sqrt(M_PI * mid));
        }

        if (sum < frac * error) {
            upper = mid;
        } else {
            lower = mid;
        }
    }

    a_vals.push_back((upper + lower) / NT(2.0));
    error = (1.0 - frac) * error;
}


// Compute a_{i+1} when a_i is given
template
<
    typename RandomPointGenerator,
    typename Polytope,
    typename Parameters,
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
NT get_next_gaussian2(Polytope &P,
                      Point &p,
                      NT a,
                      const unsigned int &N,
                      const NT &ratio,
                      const NT &C,
                      Parameters const& var,
                      RandomNumberGenerator& rng)
{

    NT last_a = a;
    NT last_ratio = 0.1;
    //k is needed for the computation of the next variance a_{i+1} = a_i * (1-1/d)^k
    NT k = 1.0;
    const NT tol = 0.00001;
    bool done=false;
    std::vector<NT> fn(N,NT(0.0));
    std::list<Point> randPoints;
    typedef typename std::vector<NT>::iterator viterator;

    //sample N points using hit and run or ball walk
    //rand_gaussian_point_generator(P, p, N, var.walk_steps, randPoints, last_a, var);

    PushBackWalkPolicy push_back_policy;
    RandomPointGenerator::apply(P, p, last_a, N, var.walk_steps, randPoints,
                                push_back_policy, rng);

    viterator fnit;
    while (!done)
    {
        a = last_a*std::pow(ratio,k);

        fnit = fn.begin();
        for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit, fnit++)
        {
            *fnit = eval_exp(*pit,a)/eval_exp(*pit, last_a);
        }
        std::pair<NT, NT> mv = getMeanVariance(fn);

        // Compute a_{i+1}
        if (mv.second/(mv.first * mv.first)>=C || mv.first/last_ratio<1.0+tol)
        {
            if (k != 1.0)
            {
                k = k / 2;
            }
            done=true;
        } else {
            k = 2 * k;
        }
        last_ratio = mv.first;
    }
    return last_a * std::pow(ratio, k);
}

// Compute the sequence of spherical gaussians
template
<
    typename RandomPointGenerator,
    typename Polytope,
    typename Parameters,
    typename NT,
    typename RandomNumberGenerator,
    typename WalkType
>
void get_annealing_schedule2(Polytope &P,
                            const NT &ratio,
                            const NT &C,
                            const NT &frac,
                            const unsigned int &N,
                            Parameters &var,
                            NT &error,
                            std::vector<NT> &a_vals,
                            RandomNumberGenerator& rng,
                            WalkType& walk)
{

    typedef typename Polytope::PointType Point;

    // Compute the first gaussian
    get_first_gaussian(P, frac, var, error, a_vals);

#ifdef VOLESTI_DEBUG
    std::cout<<"first gaussian computed\n"<<std::endl;
#endif

    NT a_stop = 0.0;
    NT curr_fn = 2.0;
    NT curr_its = 1.0;
    NT next_a;
    const NT tol = 0.001;
    unsigned int it = 0;
    unsigned int n = var.n;
    unsigned int steps;
    unsigned int coord_prev;
    const unsigned int totalSteps= ((int)150/error)+1;

    if (a_vals[0]<a_stop) a_vals[0] = a_stop;


#ifdef VOLESTI_DEBUG
    std::cout<<"Computing the sequence of gaussians..\n"<<std::endl;
#endif

    Point p(n);
    Point p_prev=p;

    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    VT lamdas;
    lamdas.setZero(P.num_of_hyperplanes());

    while (true)
    {
        if (var.ball_walk) {
            var.delta = 4.0 * var.che_rad
                    / std::sqrt(std::max(NT(1.0), a_vals[it]) * NT(n));
        }
        // Compute the next gaussian
        //next_a = get_next_gaussian2(P, p, a_vals[it], N, ratio, C, var);
        next_a = get_next_gaussian2<RandomPointGenerator>
                      (P, p, a_vals[it], N, ratio, C, var, rng);

        curr_fn = 0;
        curr_its = 0;
        lamdas.setConstant(NT(0));
        steps = totalSteps;

        /*
        if (var.cdhr_walk){
            gaussian_first_coord_point(P, p, p_prev, coord_prev, var.walk_steps,
                                       a_vals[it], lamdas, var);
            curr_its += 1.0;
            curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
            steps--;
        }*/

        // Compute some ratios to decide if this is the last gaussian
        for (unsigned  int j = 0; j < steps; j++)
        {
            //gaussian_next_point(P, p, p_prev, coord_prev, var.walk_steps,
            //                    a_vals[it], lamdas, var);

            walk.template apply(P, p, a_vals[it], var.walk_steps, rng);

            curr_its += 1.0;
            curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
        }

        // Remove the last gaussian.
        // Set the last a_i equal to zero
        if (next_a>0 && curr_fn/curr_its>(1.0+tol))
        {
            a_vals.push_back(next_a);
            it ++;
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

}

template
<
    typename RNGType,
    typename WalkTypePolicy = GaussianBallWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, double>,
    typename Polytope,
    typename GParameters,
    typename Point,
    typename NT
>
NT volume_gaussian_annealing(Polytope &P,
                             GParameters & var,
                             std::pair<Point,NT> InnerBall)
{

    typedef typename Polytope::VT 	VT;
    const NT maxNT = 1.79769e+308;
    const NT minNT = -1.79769e+308;
    NT vol;
    bool done;
    unsigned int n = var.n;
    unsigned int m = P.num_of_hyperplanes();
    unsigned int min_index, max_index, index, min_steps;
    NT error = var.error, curr_eps, min_val, max_val, val;
    NT frac = var.frac;
    typedef typename std::vector<NT>::iterator viterator;

    typedef typename WalkTypePolicy::template Walk
                                              <
                                                Polytope,
                                                RandomNumberGenerator
                                              > WalkType;
    typedef GaussianRandomPointGenerator<WalkType> RandomPointGenerator;

    RandomNumberGenerator rng(P.dimension());

    // Consider Chebychev center as an internal point
    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Save the radius of the Chebychev ball
    var.che_rad = radius;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(c.getCoefficients());

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = var.ratio;
    NT C = var.C;
    unsigned int N = var.N;

    // Computing the sequence of gaussians
#ifdef VOLESTI_DEBUG
    std::cout<<"\n\nComputing annealing...\n"<<std::endl;
    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
#endif

    WalkType walk(P, c, rng);
    get_annealing_schedule2<RandomPointGenerator>(P, ratio, C, frac,
                                                  N, var, error, a_vals, rng, walk);

#ifdef VOLESTI_DEBUG
    std::cout<<"All the variances of schedule_annealing computed in = "
            << (double)clock()/(double)CLOCKS_PER_SEC-tstart2<<" sec"<<std::endl;
    for (viterator avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++){
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;
#endif
    unsigned int mm = a_vals.size()-1;

    // Initialization for the approximation of the ratios
    unsigned int W = var.W;
    unsigned int coord_prev, i=0;
    std::vector<NT> last_W2(W,0), fn(mm,0), its(mm,0);
    VT lamdas;
    lamdas.setZero(m);
    vol=std::pow(M_PI/a_vals[0], (NT(n))/2.0);
    Point p(n), p_prev(n); // The origin is the Chebychev center of the Polytope
    viterator fnIt = fn.begin(), itsIt = its.begin(), avalsIt = a_vals.begin(), minmaxIt;

#ifdef VOLESTI_DEBUG
    std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    std::cout<<"computing ratios..\n"<<std::endl;
#endif

    // Compute the first point if CDHR is requested
    //if(var.cdhr_walk){
    //    gaussian_first_coord_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);
    //}

    for ( ; fnIt != fn.end(); fnIt++, itsIt++, avalsIt++, i++)
    { //iterate over the number of ratios
        //initialize convergence test
        curr_eps = error/std::sqrt((NT(mm)));
        done=false;
        min_val = minNT;
        max_val = maxNT;
        min_index = W-1;
        max_index = W-1;
        index = 0;
        min_steps=0;
        std::vector<NT> last_W=last_W2;

        // Set the radius for the ball walk if it is requested
        if (var.ball_walk) {
            var.delta = 4.0 * radius / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n));
        }

        while(!done || (*itsIt)<min_steps)
        {

            //gaussian_next_point(P,p,p_prev,coord_prev,var.walk_steps,*avalsIt,lamdas,var);

            walk.template apply(P, p, *avalsIt, var.walk_steps, rng);

            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
            val = (*fnIt) / (*itsIt);

            last_W[index] = val;
            if(val<=min_val){
                min_val = val;
                min_index = index;
            }else if(min_index==index){
                minmaxIt = std::min_element(last_W.begin(), last_W.end());
                min_val = *minmaxIt;
                min_index = std::distance(last_W.begin(), minmaxIt);
            }

            if(val>=max_val){
                max_val = val;
                max_index = index;
            }else if(max_index==index){
                minmaxIt = std::max_element(last_W.begin(), last_W.end());
                max_val = *minmaxIt;
                max_index = std::distance(last_W.begin(), minmaxIt);
            }

            if( (max_val-min_val)/max_val<=curr_eps/2.0 ){
                done=true;
            }

            index = index%W+1;

            if(index==W) index=0;
        }
#ifdef VOLESTI_DEBUG
        std::cout << "ratio " << i << " = " << (*fnIt) / (*itsIt)
                  << " N_" << i << " = " << *itsIt << std::endl;
#endif
        vol = vol*((*fnIt) / (*itsIt));
    }

#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
#endif

    P.free_them_all();
    return vol;
}


#endif
