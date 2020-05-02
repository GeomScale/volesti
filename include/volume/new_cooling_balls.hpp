// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

#ifndef NEW_COOLING_BALLS_H
#define NEW_COOLING_BALLS_H

#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "ballintersectconvex.h"
#include "vpolyintersectvpoly.h"
//#include "samplers.h"
#include "rounding.h"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "ball_annealing.h"
#include "ratio_estimation.h"


////////////////////////////////////
// ball annealing

template <typename NT>
struct cooling_ball_parameters
{
    cooling_ball_parameters()
        :   lb(0.1)
        ,   ub(0.15)
        ,   p(0.75)
        ,   rmax(0)
        ,   alpha(0.2)
        ,   win_len(500)
        ,   N(150)
        ,   nu(10)
        ,   window2(false)
    {}

    NT lb;
    NT ub;
    NT p;
    NT rmax;
    NT alpha;
    int win_len;
    int N;
    int nu;
    bool window2;
};

/// Helpers

template <typename Point, typename ConvexBody, typename PointList, typename NT>
bool check_convergence(ConvexBody const& P,
                       PointList const& randPoints,
                       bool& too_few,
                       NT& ratio,
                       int const& nu,
                       bool const& precheck,
                       bool const& lastball,
                       cooling_ball_parameters<NT> const& parameters)
{
    NT alpha = parameters.alpha;
    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int m = randPoints.size()/nu;
    int i = 1;
    NT T;
    NT rs;
    NT alpha_check = 0.01;
    size_t countsIn = 0;

    for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++)
    {
        if (P.is_in(*pit)==-1) countsIn++;
        if (i % m == 0)
        {
            ratios.push_back(NT(countsIn)/m);
            countsIn = 0;
            if (ratios.size()>1 && precheck)
            {
                boost::math::students_t dist(ratios.size() - 1);
                mv = get_mean_variance(ratios);
                ratio = mv.first;
                rs = std::sqrt(mv.second);
                T = rs * (boost::math::quantile
                            (boost::math::complement(dist, alpha_check / 2.0))
                          / std::sqrt(NT(ratios.size())));
                if (ratio + T < parameters.lb)
                {
                    too_few = true;
                    return false;
                } else if (ratio - T > parameters.ub) return false;
            }
        }
    }

    if (precheck) alpha *= 0.5;
    mv = get_mean_variance(ratios);
    ratio = mv.first;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs * (boost::math::quantile(boost::math::complement(dist, alpha))
           / std::sqrt(NT(nu)));
    if (ratio > parameters.lb + T)
    {
        if (lastball) return true;
        if ((precheck && ratio < parameters.ub - T)
        || (!precheck && ratio < parameters.ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;
}


template <typename Polytope, typename Ball, typename NT, typename RNG>
bool get_first_ball(Polytope const& P,
                    Ball& B0,
                    NT& ratio,
                    NT const& radius_input,
                    cooling_ball_parameters<NT> const& parameters,
                    RNG& rng)
{
    const unsigned max_iterarions = 10;
    const unsigned tolerance = 0.00000000001;
    typedef typename Polytope::PointType Point;
    int n = P.dimension();
    int iter = 1;
    bool bisection_int = false;
    bool pass = false;
    bool too_few = false;
    std::list<Point> randPoints;
    NT rmax = parameters.rmax;
    NT sqrt_n = std::sqrt(NT(n));
    NT rad1 = radius_input;

    if (rmax > 0.0)
    {
        for (int i = 0; i < 1200; ++i)
        {
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rmax, rng));
        }
        pass = check_convergence<Point>(P, randPoints, too_few, ratio,
                                        10, true, false, parameters);
        if (pass || !too_few)
        {
            B0 = Ball(Point(n), rmax * rmax);
            return true;
        }
        bisection_int = true;
    } else
    {
        rmax = 2 * sqrt_n * rad1;
    }
    NT radius = rad1;

    while (!bisection_int)
    {
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i)
        {
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rmax, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters))
        {
            B0 = Ball(Point(n), rmax * rmax);
            return true;
        }

        if (too_few) break;
        rad1 = rmax;
        rmax = rmax + 2 * sqrt_n * radius;
    }

    NT rad_med;
    NT rad0=rad1;
    NT rad_m = rmax;

    while (iter <= max_iterarions)
    {
        rad_med = 0.5*(rad1+rmax);
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i)
        {
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rad_med, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters))
        {
            B0 = Ball(Point(n), rad_med * rad_med);
            return true;
        }

        if (too_few)
        {
            rmax = rad_med;
        } else {
            rad1 = rad_med;
        }

        if (rmax-rad1 < tolerance)
        {
            rad1 = rad0;
            rmax = rad_m;
            iter++;
        }
    }
    return false;
}


template <typename Point, typename ball, typename PointList, typename NT>
bool get_next_zonotopeball(std::vector<ball>& BallSet,
                           PointList const& randPoints,
                           NT const& rad_min,
                           std::vector<NT>& ratios,
                           cooling_ball_parameters<NT> const& parameters)
{
    const unsigned max_iterarions = 10;
    const unsigned tolerance = 0.00000000001;
    int n = (*randPoints.begin()).dimension();
    int iter = 1;
    bool too_few;
    NT radmax = NT(0);
    NT radmin = rad_min;

    for (auto rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit)
    {
        NT pnorm = (*rpit).squared_length();
        if (pnorm > radmax) radmax = pnorm;
    }
    ball Biter;
    radmax = std::sqrt(radmax);
    NT radmin_init = radmin;
    NT radmax_init = radmax;

    while (iter <= max_iterarions)
    {
        NT rad = 0.5 * (radmin + radmax);
        Biter = ball(Point(n), rad * rad);
        too_few = false;

        NT ratio;
        if (check_convergence<Point>(Biter, randPoints, too_few, ratio,
                                     parameters.nu, false, false, parameters))
        {
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return true;
        }

        if (too_few)
        {
            radmin = rad;
        } else
        {
            radmax = rad;
        }

        if (radmax-radmin < tolerance)
        {
            radmin = radmin_init;
            radmax = radmax_init;
            iter++;
        }
    }
    return false;
}


template
<
    typename RandomPointGenerator,
    typename PolyBall,
    typename ball,
    typename Polytope,
    typename NT,
    typename RNG
>
bool get_sequence_of_polytopeballs(Polytope const& P,
                                   std::vector<ball>& BallSet,
                                   std::vector<NT>& ratios,
                                   int const& Ntot,
                                   NT const& radius,
                                   unsigned int const& walk_length,
                                   NT& diameter,
                                   cooling_ball_parameters<NT> parameters,
                                   RNG& rng)
{

    typedef typename Polytope::PointType Point;
    bool fail;
    int n = P.dimension();
    NT ratio;
    NT ratio0;
    std::list<Point> randPoints;
    ball B0;
    Point q(n);
    PolyBall zb_it;

    if ( !get_first_ball(P, B0, ratio, radius, parameters, rng) )
    {
        return false;
    }

    ratio0 = ratio;

    PushBackWalkPolicy push_back_policy;
    RandomPointGenerator::apply(P, q, Ntot, walk_length,
                                randPoints, push_back_policy, rng);

    if (check_convergence<Point>(B0, randPoints,
                                 fail, ratio, parameters.nu,
                                 false, true, parameters))
    {
        ratios.push_back(ratio);
        BallSet.push_back(B0);
        ratios.push_back(ratio0);
        return true;
    }
    if ( !get_next_zonotopeball<Point>(BallSet, randPoints, B0.radius(), ratios,
                                       parameters) )
    {
        return false;
    }

    while (true)
    {
        zb_it = PolyBall(P, BallSet[BallSet.size()-1]);
        q = Point(n);
        randPoints.clear();
        zb_it.comp_diam(diameter);

        RandomPointGenerator::apply(zb_it, q, Ntot, walk_length,
                                    randPoints, push_back_policy, rng);
        if (check_convergence<Point>(B0, randPoints, fail, ratio, parameters.nu,
                                     false, true, parameters))
        {
            ratios.push_back(ratio);
            BallSet.push_back(B0);
            ratios.push_back(ratio0);
            return true;
        }
        if ( !get_next_zonotopeball<Point>(BallSet, randPoints, B0.radius(),
                                           ratios, parameters) )
        {
            return false;
        }
    }
}


////////////////////////////////////
///
/// ratio estimation

template <typename NT>
bool is_max_error(NT const& a, NT const& b, NT const& error)
{
    return ((b-a)/a<error/2.0) ? true : false;
}

template
<
    typename WalkType,
    typename Point,
    typename PolyBall1,
    typename PolyBall2,
    typename NT,
    typename RNG
>
NT estimate_ratio(PolyBall1 const& Pb1,
                  PolyBall2 const& Pb2,
                  NT const& ratio,
                  NT const& error,
                  int const& W,
                  int const& Ntot,
                  unsigned int const& walk_length,
                  RNG& rng,
                  bool isball = false,
                  NT radius = 0.0)
{
    const unsigned max_iterations_estimation = 10000000;
    int n = Pb1.dimension();
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int iter = 1;
    NT min_val = std::numeric_limits<NT>::lowest();
    NT max_val = std::numeric_limits<NT>::max();
    NT val;
    size_t totCount = Ntot;
    size_t countIn = Ntot * ratio;
    std::vector<NT> last_W(W);

    typename std::vector<NT>::iterator minmaxIt;
    Point p(n);
    WalkType walk(Pb1, p, rng);

    while (iter++ <= max_iterations_estimation)
    {
        if (isball)
        {
            p = GetPointInDsphere<Point>::apply(n, radius, rng);
        } else {
            walk.template apply(Pb1, p, walk_length, rng);
        }
        if (Pb2.is_in(p)==-1) countIn = countIn + 1.0;

        totCount = totCount + 1.0;
        val = NT(countIn) / NT(totCount);
        last_W[index] = val;

        if (val<=min_val)
        {
            min_val = val;
            min_index = index;
        } else if (min_index==index)
        {
            minmaxIt = std::min_element(last_W.begin(), last_W.end());
            min_val = *minmaxIt;
            min_index = std::distance(last_W.begin(), minmaxIt);
        }

        if (val>=max_val)
        {
            max_val = val;
            max_index = index;
        } else if (max_index==index)
        {
            minmaxIt = std::max_element(last_W.begin(), last_W.end());
            max_val = *minmaxIt;
            max_index = std::distance(last_W.begin(), minmaxIt);
        }

        if ( (max_val-min_val)/max_val<=error/2.0 )
        {
            return val;
        }

        index = index%W + 1;
        if (index==W) index=0;
    }
    return val;
}


template
<
    typename WalkType,
    typename Point,
    typename PolyBall1,
    typename PolyBall2,
    typename NT,
    typename RNG
>
NT estimate_ratio_interval(PolyBall1 const& Pb1,
                           PolyBall2 const& Pb2,
                           NT const& ratio,
                           NT const& error,
                           int const& W,
                           int const& Ntot,
                           NT const& prob,
                           unsigned int const& walk_length,
                           RNG& rng,
                           bool isball = false,
                           NT const& radius = 0)
{
    const unsigned max_iterations_estimation = 10000000;
    int n = Pb1.dimension();
    int index = 0;
    int iter = 1;
    std::vector<NT> last_W(W);
    NT val;
    NT sum_sq = NT(0);
    NT sum = NT(0);
    size_t totCount = Ntot;
    size_t countIn = Ntot * ratio;
    Point p(n);
    WalkType walk(Pb1, p, rng);

    for (int i = 0; i < W; ++i)
    {
        if (isball)
        {
            p = GetPointInDsphere<Point>::apply(n, radius, rng);
        } else
        {
            walk.template apply(Pb1, p, walk_length, rng);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1;

        totCount = totCount + 1;
        val = NT(countIn) / NT(totCount);
        sum += val;
        sum_sq += val * val;
        last_W[index] = val;
        index = index % W + 1;
        if (index == W) index = 0;
    }

    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0));
    NT m=sum/NT(W);
    NT s;

    while (iter++ <= max_iterations_estimation)
    {
        if (isball) {
            p = GetPointInDsphere<Point>::apply(n, radius, rng);
        } else
        {
            walk.template apply(Pb1, p, walk_length, rng);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1;

        totCount = totCount + 1;
        val = NT(countIn) / NT(totCount);

        m = (m - last_W[index] / NT(W)) + val / NT(W);
        sum_sq = (sum_sq - last_W[index] * last_W[index]) + val * val;
        sum = (sum - last_W[index]) + val;
        s = std::sqrt((sum_sq + NT(W) * m * m - 2.0 * m * sum) / NT(W));
        last_W[index] = val;

        index = index % W + 1;
        if (index == W) index = 0;

        if (is_max_error(val - zp * s, val + zp * s, error))
        {
            return val;
        }
    }
    return val;
}



template
<
    typename WalkTypePolicy = BilliardWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, double>,
    typename Polytope
>
double volume_cooling_balls(Polytope const& Pin,
                            double const& error = 1.0,
                            unsigned int const& walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope <Polytope, BallType> PolyBall;
    typedef typename Polytope::VT VT;
    typedef std::list <Point> PointList;

    typedef typename WalkTypePolicy::template Walk
                                              <
                                                Polytope,
                                                RandomNumberGenerator
                                              > WalkType;
    typedef RandomPointGenerator<WalkType> RandomPointGenerator;


    auto P(Pin);
    RandomNumberGenerator rng(P.dimension());
    cooling_ball_parameters<NT> parameters;

    int n = P.dimension();
    NT prob = parameters.p;
    int N_times_nu = parameters.N * parameters.nu;


    auto InnerBall = P.ComputeInnerBall();
    NT radius = InnerBall.second;
    Point c = InnerBall.first;
    NT diameter = P.ComputeDiameter();

    std::vector<BallType> BallSet;
    std::vector<NT> ratios;

    // Normalize and move the chebychev center to the origin
    // and apply the same shifting to the polytope
    P.normalize();
    P.shift(c.getCoefficients());

    //Point new_c(n); //origin
    //WalkType walk(P, new_c, rng);

    if ( !get_sequence_of_polytopeballs
          <
            RandomPointGenerator,
            PolyBall
          >(P, BallSet, ratios,
            N_times_nu, radius, walk_length, diameter,
            parameters, rng) )
    {
        return -1.0;
    }

    NT vol = (std::pow(M_PI, n / 2.0)
              * (std::pow((*(BallSet.end() - 1)).radius(), n)))
            / (tgamma(n / 2.0 + 1));

    int mm = BallSet.size() + 1;
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = error / (2.0 * std::sqrt(NT(mm)));
    NT er1 = (error * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    vol *= (parameters.window2) ?
                estimate_ratio<WalkType, Point>(*(BallSet.end() - 1),
                                      P,
                                      *(ratios.end() - 1),
                                      er0, parameters.win_len, 1200, walk_length,
                                      rng,
                                      true,
                                      (*(BallSet.end() - 1)).radius())
              : estimate_ratio_interval<WalkType, Point>(*(BallSet.end() - 1),
                                               P,
                                               *(ratios.end() - 1),
                                               er0, parameters.win_len, 1200, prob,
                                               walk_length, rng,
                                               true,
                                               (*(BallSet.end() - 1)).radius());

    PolyBall Pb;
    auto balliter = BallSet.begin();
    auto ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1)
    {
        vol *= (!parameters.window2) ?
               1 / estimate_ratio_interval<WalkType, Point>(P, *balliter, *ratioiter,
                                                  er1,
                                                  parameters.win_len,
                                                  N_times_nu,
                                                  prob,
                                                  walk_length, rng)
            : 1 / estimate_ratio<WalkType, Point>(P, *balliter, *ratioiter,
                                        er1, parameters.win_len,
                                        N_times_nu,
                                        walk_length, rng);
    }
    for ( ; balliter < BallSet.end() - 1; ++balliter, ++ratioiter)
    {
        Pb = PolyBall(P, *balliter);
        Pb.comp_diam(diameter);
        vol *= (!parameters.window2) ?
                    1 / estimate_ratio_interval<WalkType, Point>(Pb,
                                                       *(balliter + 1),
                                                       *(ratioiter + 1),
                                                       er1, parameters.win_len,
                                                       N_times_nu,
                                                       prob, walk_length,
                                                       rng)
                  : 1 / estimate_ratio<WalkType, Point>(Pb, *balliter, *ratioiter, er1,
                                              parameters.win_len,
                                              N_times_nu,
                                              walk_length, rng);
    }

    P.free_them_all();
    return vol;
}

#endif
