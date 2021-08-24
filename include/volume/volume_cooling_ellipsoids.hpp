// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar

// Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_COOLING_ELLIPSOIDS_HPP
#define VOLUME_COOLING_ELLIPSOIDS_HPP

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/orderpolytope.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/ellipsoid.h"
#include "sampling/random_point_generators.hpp"
#include "random_walks/gaussian_accelerated_billiard_walk.hpp"
#include "preprocess/max_inscribed_ellipsoid.hpp"
#include "volume/math_helpers.h"
#include "volume/sampling_policies.hpp"
#include "volume/volume_cooling_balls.hpp"


////////////////////////////////////
// ellipsoid annealing

template <typename NT>
struct cooling_ellipsoid_parameters
{
    cooling_ellipsoid_parameters(unsigned int win_len)
        :   lb(0.1)
        ,   ub(0.15)
        ,   p(0.75)
        ,   q_max(0)
        ,   alpha(0.2)
        ,   win_len(win_len)
        ,   N(125)
        ,   nu(10)
        ,   window2(false)
    {}

    NT lb;
    NT ub;
    NT p;
    NT q_max;
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
                       bool const& last_ellipsoid,
                       cooling_ellipsoid_parameters<NT> const& parameters)
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
        if (last_ellipsoid) return true;
        if ((precheck && ratio < parameters.ub - T)
        || (!precheck && ratio < parameters.ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;
}


template <typename Polytope, typename Ellipsoid, typename NT, typename RNG>
std::pair<NT, NT> get_first_ellipsoid(Polytope const& P,
                    Ellipsoid& E0,
                    NT& ratio,
                    Ellipsoid const& inscribed_ellipsoid,
                    cooling_ellipsoid_parameters<NT> const& parameters,
                    RNG& rng) {
    typedef typename Ellipsoid::VT VT;
    typedef typename Ellipsoid::MT MT;

    const unsigned max_iterarions = 20;
    NT tolerance = 0.00000000001;
    typedef typename Polytope::PointType Point;
    int n = P.dimension();
    int iter = 1;
    bool bisection_int = false;
    bool pass = false;
    bool too_few = false;
    std::list <Point> randPoints;
    NT q_max = parameters.q_max;
    NT sqrt_n = std::sqrt(NT(n));
    NT q1 = 1.0;
    q_max = 2 * sqrt_n;   // 2 * sqrt_n times the inscribed ellipsoid

    std::cout << "in_ell radius: " << inscribed_ellipsoid.radius() << std::endl;
    E0 = inscribed_ellipsoid;

    while (!bisection_int) {
        randPoints.clear();
        too_few = false;

        E0.scale(q_max/q1);
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(GetPointInDellipsoid<Point>::apply(n, E0, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters)) {
            return std::pair<NT, NT> (q_max, q_max);
        }

        if (too_few) break;
        q1 = q_max;
        q_max = q_max + 2 * sqrt_n;
    }

    NT q_mid;
    NT q0 = q1;
    NT q_m = q_max;
    NT q_scale_prev = q_max;

    while (iter <= max_iterarions) {
        q_mid = 0.5 * (q1 + q_max);
        randPoints.clear();
        too_few = false;

        E0.scale(q_mid/q_scale_prev);
        q_scale_prev = q_mid;
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(GetPointInDellipsoid<Point>::apply(n, E0, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters)) {
            return std::pair<NT, NT> (q_mid, q_m);
        }

        if (too_few) {
            q_max = q_mid;
        } else {
            q1 = q_mid;
        }

        if (q_max - q1 < tolerance) {
            q1 = q0;
            q_max = q_m;
            iter++;
        }
    }

    return std::pair<NT, NT> (-1, -1);
}


template <typename Point, typename Ellipsoid, typename PointList, typename NT>
NT get_next_ellipsoid(std::vector<Ellipsoid>& EllipsoidSet,
                           PointList const& randPoints,
                           Ellipsoid const& inscribed_ellipsoid,
                           NT q_min,
                           NT q_max,
                           std::vector<NT>& ratios,
                           cooling_ellipsoid_parameters<NT> const& parameters)
{
    const unsigned max_iterarions = 20;
    NT tolerance = 0.00000000001;
    int n = (*randPoints.begin()).dimension();
    int iter = 1;
    bool too_few;

    Ellipsoid Eiter = inscribed_ellipsoid;
    std::cout << "q_min: " << q_min << ", q_max: " << q_max << std::endl;
    NT q_min_init = q_min;
    NT q_max_init = q_max;

    NT q_prev = 1.0;
    while (iter <= max_iterarions)
    {
        NT q = 0.5 * (q_min + q_max);
        Eiter.scale(q/q_prev);
        q_prev = q;
        too_few = false;

        NT ratio;
        if (check_convergence<Point>(Eiter, randPoints, too_few, ratio,
                                     parameters.nu, false, false, parameters))
        {
            EllipsoidSet.push_back(Eiter);
            ratios.push_back(ratio);
            std::cout << "selected q: " << q << std::endl;
            std::cout << "E_iter radius: " << Eiter.radius() << std::endl;
            return q;
        }

        if (too_few)
        {
            q_min = q;
        } else
        {
            q_max = q;
        }

        if (q_max-q_min < tolerance)
        {
            q_min = q_min_init;
            q_max = q_max_init;
            iter++;
        }
    }

    return NT(-1);
}


template
<
    typename RandomPointGenerator,
    typename PolyEllipsoid,
    typename Ellipsoid,
    typename Polytope,
    typename NT,
    typename RNG
>
bool get_sequence_of_polytope_ellipsoids(Polytope& P,
                                   std::vector<Ellipsoid>& EllipsoidSet,
                                   std::vector<NT>& ratios,
                                   int const& Ntot,
                                   Ellipsoid const& inscribed_ellipsoid,
                                   unsigned int const& walk_length,
                                   cooling_ellipsoid_parameters<NT> const& parameters,
                                   RNG& rng)
{

    typedef typename Polytope::PointType Point;
    bool fail;
    int n = P.dimension();
    NT ratio;
    NT ratio0;
    std::list<Point> randPoints;
    Ellipsoid E0;
    Point x0(n);

    std::pair<NT, NT> q_init = get_first_ellipsoid(P, E0, ratio, inscribed_ellipsoid, parameters, rng);
    if (q_init.first == -1)
    {
        return false;
    }

    ratio0 = ratio;
    PushBackWalkPolicy push_back_policy;
    RandomPointGenerator::apply(P, x0, inscribed_ellipsoid, Ntot, walk_length,
                                randPoints, push_back_policy, rng);

    std::cout << "radius of E0: " << E0.radius() << std::endl;
    if (check_convergence<Point>(E0, randPoints,
                                 fail, ratio, parameters.nu,
                                 false, true, parameters))
    {
        ratios.push_back(ratio);
        EllipsoidSet.push_back(E0);
        ratios.push_back(ratio0);
        return true;
    }

    // for further ellipsoids
    NT q_min = q_init.first;
    NT q_max = q_init.second;

    q_max = get_next_ellipsoid<Point>(EllipsoidSet, randPoints, inscribed_ellipsoid, q_min, q_max, ratios, parameters);
    if (q_max == -1)
    {
        return false;
    }
    std::cout << "Got second ellipsoid" << std::endl;

    while (true)
    {
        PolyEllipsoid zb_it(P, EllipsoidSet[EllipsoidSet.size()-1]);
        x0.set_to_origin();
        randPoints.clear();

        RandomPointGenerator::apply(zb_it, x0, inscribed_ellipsoid, Ntot, walk_length,
                                    randPoints, push_back_policy, rng);
        if (check_convergence<Point>(E0, randPoints, fail, ratio, parameters.nu,
                                     false, true, parameters))
        {
            ratios.push_back(ratio);
            EllipsoidSet.push_back(E0);
            ratios.push_back(ratio0);
            return true;
        }

        q_max = get_next_ellipsoid<Point>(EllipsoidSet, randPoints, inscribed_ellipsoid, q_min, q_max, ratios, parameters);
        if (q_max == -1)
        {
            return false;
        }
    }
}


template
<
        typename WalkType,
        typename Point,
        typename PolyEllipsoid1,
        typename PolyEllipsoid2,
        typename Ellipsoid,
        typename NT,
        typename RNG
>
NT estimate_ratio_ellipsoid(PolyEllipsoid1 const& Pe1,
                            PolyEllipsoid2 const& Pe2,
                            Ellipsoid const& E,
                            NT const& ratio,
                            NT const& error,
                            unsigned int const& W,
                            unsigned int const& Ntot,
                            unsigned int const& walk_length,
                            RNG& rng)
{
    estimate_ratio_parameters<NT> ratio_parameters(W, Ntot, ratio);
    unsigned int n = Pe1.dimension();
    Point p(n);
    WalkType walk(Pe1, p, E, rng);

    do
    {
        walk.template apply(Pe1, p, E, walk_length, rng);
    } while(!estimate_ratio_generic(Pe2, p, error, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}

template
<       typename Point,
        typename Ellipsoid,
        typename PolyEllipsoid,
        typename NT,
        typename RNG
>
NT estimate_ratio_ellipsoid(Ellipsoid const& E,
                            PolyEllipsoid const& Pe2,
                            NT const& ratio,
                            NT const& error,
                            int const& W,
                            int const& Ntot,
                            RNG& rng)
{
    estimate_ratio_parameters<NT> ratio_parameters(W, Ntot, ratio);
    unsigned int n = E.dimensions();
    Point p(n);

    do
    {
        p = GetPointInDellipsoid<Point>::apply(n, E, rng);
    } while(!estimate_ratio_generic(Pe2, p, error, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}


template
<
    typename Point,
    typename Ellipsoid,
    typename PolyEllipsoid2,
    typename NT,
    typename RNG
>
NT estimate_ratio_interval_ellipsoid(Ellipsoid const& E,
                                     PolyEllipsoid2 const& Pe2,
                                     NT const& ratio,
                                     NT const& error,
                                     int const& W,
                                     int const& Ntot,
                                     NT const& prob,
                                     RNG& rng)
{
    estimate_ratio_interval_parameters<NT> ratio_parameters(W, Ntot, ratio);
    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0));

    unsigned int n = Pe2.dimension();
    Point p(n);

    for (int i = 0; i < ratio_parameters.W; ++i)
    {
        p = GetPointInDellipsoid<Point>::apply(n, E, rng);
        full_sliding_window(Pe2, p, ratio_parameters);
    }

    ratio_parameters.mean = ratio_parameters.sum / NT(ratio_parameters.W);

    do {
        p = GetPointInDellipsoid<Point>::apply(n, E, rng);
    } while (!estimate_ratio_interval_generic(Pe2, p, error, zp, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}


template
<
        typename WalkType,
        typename Point,
        typename PolyEllipsoid1,
        typename PolyEllipsoid2,
        typename Ellipsoid,
        typename NT,
        typename RNG
>
NT estimate_ratio_interval_ellipsoid(PolyEllipsoid1 const& Pe1,
                                     PolyEllipsoid2 const& Pe2,
                                     Ellipsoid const& E,
                                     NT const& ratio,
                                     NT const& error,
                                     int const& W,
                                     int const& Ntot,
                                     NT const& prob,
                                     unsigned int const& walk_length,
                                     RNG& rng)
{
    estimate_ratio_interval_parameters<NT> ratio_parameters(W, Ntot, ratio);
    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0));

    unsigned int n = Pe1.dimension();
    Point p(n);
    WalkType walk(Pe1, p, E, rng);

    for (int i = 0; i < ratio_parameters.W; ++i)
    {
        walk.template apply(Pe1, p, E, walk_length, rng);
        full_sliding_window(Pe2, p, ratio_parameters);
    }

    ratio_parameters.mean = ratio_parameters.sum / NT(ratio_parameters.W);

    do {
        walk.template apply(Pe1, p, E, walk_length, rng);
    } while (!estimate_ratio_interval_generic(Pe2, p, error, zp, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}


template
<
    typename WalkTypePolicy,
    typename Polytope,
    typename RandomNumberGenerator
>
std::pair<double, double> volume_cooling_ellipsoids(Polytope const& Pin,
                                       RandomNumberGenerator &rng,
                                       double const& error = 0.1,
                                       unsigned int const& walk_length = 1,
                                       unsigned int const& win_len = 600)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ellipsoid<Point> EllipsoidType;
    typedef BallIntersectPolytope <Polytope, EllipsoidType> PolyEllipsoid;
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    typedef std::list <Point> PointList;

    typedef typename WalkTypePolicy::template Walk
                                              <
                                                Polytope,
                                                RandomNumberGenerator
                                              > WalkType;
    typedef MultivariateGaussianRandomPointGenerator<WalkType> RandomPointGenerator;

    auto P(Pin);
    cooling_ellipsoid_parameters<NT> parameters(win_len);

    int n = P.dimension();
    NT prob = parameters.p;
    int N_times_nu = parameters.N * parameters.nu;

    // ----------- Get inscribed ellipsoid --------------------------------
    unsigned int max_iter = 150;
    NT tol = std::pow(10, -6.0), reg = std::pow(10, -4.0);

    std::pair<Point, NT> InnerBall = P.ComputeInnerBall();
    if (InnerBall.second < 0.0) return std::pair<NT, NT> (-1.0, 0.0);
    VT x0 = InnerBall.first.getCoefficients();
    std::pair<std::pair<MT, VT>, bool> inscribed_ellipsoid_res = max_inscribed_ellipsoid<MT>(P.get_mat(),
                                                                                             P.get_vec(),
                                                                                             x0,
                                                                                             max_iter,
                                                                                             tol,
                                                                                             reg);
    if (!inscribed_ellipsoid_res.second)        // not converged
        return std::pair<NT, NT> (-1.0, 0.0);

    MT A_ell = inscribed_ellipsoid_res.first.first.inverse();
    EllipsoidType inscribed_ellipsoid( A_ell );
    MT L_cov = inscribed_ellipsoid.Lcov();
    Point c = inscribed_ellipsoid_res.first.second;
    // --------------------------------------------------------------------
    // // ---------------------- Inner ball as inscribed_ellipsoid --------------
    // std::pair<Point, NT> InnerBall = P.ComputeInnerBall();
    // if (InnerBall.second < 0.0) return std::pair<NT, NT> (-1.0, 0.0);
    // Point c = InnerBall.first;
    // MT A_ell = (NT(1.0)/(InnerBall.second*InnerBall.second)) * Eigen::MatrixXd::Identity(P.dimension(), P.dimension());
    // EllipsoidType inscribed_ellipsoid(A_ell);
    // // -----------------------------------------------------------------------

    std::vector<EllipsoidType> EllipsoidSet;
    std::vector<NT> ratios;

    // move the chebychev center to the origin
    // and apply the same shifting to the polytope
    P.shift(c.getCoefficients());

    if ( !get_sequence_of_polytope_ellipsoids
          <
            RandomPointGenerator,
            PolyEllipsoid
          >(P, EllipsoidSet, ratios,
            N_times_nu, inscribed_ellipsoid, walk_length,
            parameters, rng) )
    {
        return std::pair<NT, NT> (-1.0, 0.0);
    }

    std::cout << "Number of ellipsoids: " << EllipsoidSet.size() << std::endl;


    NT vol = (*(EllipsoidSet.end() - 1)).log_volume();
    std::cout << "vol: " << vol << std::endl;

    int mm = EllipsoidSet.size() + 1;
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = error / (2.0 * std::sqrt(NT(mm)));
    NT er1 = (error * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    vol += (parameters.window2) ?
                std::log(estimate_ratio_ellipsoid<Point>(*(EllipsoidSet.end() - 1),
                                                         P, *(ratios.end() - 1),
                                                         er0, parameters.win_len, 1200, rng))
              : std::log(estimate_ratio_interval_ellipsoid<Point>(*(EllipsoidSet.end() - 1),
                                                                  P, *(ratios.end() - 1),
                                                                  er0, parameters.win_len, 1200,
                                                                  prob, rng));
    std::cout << "vol: " << vol << std::endl;

    auto ellipsoiditer = EllipsoidSet.begin();
    auto ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1)
    {
        std::cout << "reached ratioiter != 1 " << std::endl;
        vol += (!parameters.window2) ?
               std::log(NT(1) / estimate_ratio_interval_ellipsoid
                    <WalkType, Point>(P,
                                      *ellipsoiditer,
                                      inscribed_ellipsoid,
                                      *ratioiter,
                                      er1,
                                      parameters.win_len,
                                      N_times_nu,
                                      prob,
                                      walk_length,
                                      rng))
            : std::log(NT(1) / estimate_ratio_ellipsoid
                    <WalkType, Point>(P,
                                      *ellipsoiditer,
                                      inscribed_ellipsoid,
                                      *ratioiter,
                                      er1,
                                      parameters.win_len,
                                      N_times_nu,
                                      walk_length,
                                      rng));
    }
    std::cout << "vol: " << vol << std::endl;

    for ( ; ellipsoiditer < EllipsoidSet.end() - 1; ++ellipsoiditer, ++ratioiter)
    {
        PolyEllipsoid Pe(P, *ellipsoiditer);
        vol += (!parameters.window2) ?
                    std::log(NT(1) / estimate_ratio_interval_ellipsoid
                                <WalkType, Point>(Pe,
                                                  *(ellipsoiditer + 1),
                                                  inscribed_ellipsoid,
                                                  *(ratioiter + 1),
                                                  er1, parameters.win_len,
                                                  N_times_nu,
                                                  prob, walk_length,
                                                  rng))
                  : std::log(NT(1) / estimate_ratio_ellipsoid
                                <WalkType, Point>(Pe,
                                                  *ellipsoiditer,
                                                  inscribed_ellipsoid,
                                                  *ratioiter,
                                                  er1,
                                                  parameters.win_len,
                                                  N_times_nu,
                                                  walk_length,
                                                  rng));
        std::cout << "vol: " << vol << std::endl;
    }

    return std::pair<NT, NT> (vol, std::exp(vol));
}



template
<
    typename WalkTypePolicy = GaussianAcceleratedBilliardWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt11213b,
                                                                double>,
    typename Polytope
>
std::pair<double, double> volume_cooling_ellipsoids(Polytope const& Pin,
                                               double const& error = 0.1,
                                               unsigned int const& walk_length = 1,
                                               unsigned int const& win_length = 600)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_cooling_ellipsoids<WalkTypePolicy>(Pin, rng, error, walk_length, win_length);
}


#endif // VOLUME_COOLING_ELLIPSOIDS_HPP