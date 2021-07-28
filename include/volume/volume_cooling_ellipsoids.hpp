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

#include "random_walks.hpp"
#include "preprocess/max_inscribed_ellipsoid.hpp"
#include "math_helpers.h"
#include "volume_cooling_balls.hpp"


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
                    // ! NT const& radius_input,
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
    // ! NT radius1 = radius_input;

    VT eigenvals_inv_sqrt = inscribed_ellipsoid.eigenvals_inv_sqrt();
    MT eigenvecs = inscribed_ellipsoid.eigenvecs();

    q_max = 2 * sqrt_n;   // 2 * sqrt_n times the inscribed ellipsoid
    VT temp_eigenvals_inv_sqrt = eigenvals_inv_sqrt;

    E0 = inscribed_ellipsoid;

    while (!bisection_int) {
        randPoints.clear();
        too_few = false;

        temp_eigenvals_inv_sqrt = q_max * eigenvals_inv_sqrt;
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(GetPointInDellipsoid<Point>::apply(n, temp_eigenvals_inv_sqrt, eigenvecs, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters)) {
            // ! B0 = Ball(Point(n), q_max * q_max);
            E0.scale(q_max);
            return std::pair<NT, NT> (q_max, q_max);
            // ! return true;
        }

        if (too_few) break;
        q1 = q_max;
        q_max = q_max + 2 * sqrt_n;
    }

    NT q_mid;
    NT q0 = q1;
    NT q_m = q_max;

    while (iter <= max_iterarions) {
        q_mid = 0.5 * (q1 + q_max);
        randPoints.clear();
        too_few = false;

        temp_eigenvals_inv_sqrt = q_mid * eigenvals_inv_sqrt;
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(GetPointInDellipsoid<Point>::apply(n, temp_eigenvals_inv_sqrt, eigenvecs, rng));
        }

        if (check_convergence<Point>(P, randPoints, too_few, ratio, 10,
                                     true, false, parameters)) {
            // ! B0 = Ball(Point(n), q_mid * q_mid);
            E0.scale(q_mid);
            return std::pair<NT, NT> (q_mid, q_max);
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
                        // !   NT const& rad_min,
                           std::vector<NT>& ratios,
                           cooling_ellipsoid_parameters<NT> const& parameters)
{
    const unsigned max_iterarions = 20;
    NT tolerance = 0.00000000001;
    int n = (*randPoints.begin()).dimension();
    int iter = 1;
    bool too_few;
    // ! NT q_max = NT(0);
    // ! NT radmin = rad_min;
    // ! NT q_min = 1.0;

    // ! for (auto rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit)
    // ! {
    // !    NT pnorm = (*rpit).squared_length();
    // !    if (pnorm > radmax) q_max = pnorm;
    // ! }

    Ellipsoid Eiter = inscribed_ellipsoid;
    // ! q_max = std::sqrt(q_max);
    NT q_min_init = q_min;
    NT q_max_init = q_max;

    NT q_prev = 1.0;
    while (iter <= max_iterarions)
    {
        NT q = 0.5 * (q_min + q_max);
        Eiter.scale(q/q_prev);
        q_prev = q;

        // ! Eiter = ellipsoid(Point(n), q * q);     // TODO: Create new ellipsoid
        too_few = false;

        NT ratio;
        if (check_convergence<Point>(Eiter, randPoints, too_few, ratio,
                                     parameters.nu, false, false, parameters))
        {
            EllipsoidSet.push_back(Eiter);
            ratios.push_back(ratio);
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

    typename Ellipsoid::MT L_cov = inscribed_ellipsoid.Lcov();
    PushBackWalkPolicy push_back_policy;
    RandomPointGenerator::apply(P, x0, L_cov, Ntot, walk_length,
                                randPoints, push_back_policy, rng);

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

    while (true)
    {
        PolyEllipsoid zb_it(P, EllipsoidSet[EllipsoidSet.size()-1]);
        x0.set_to_origin();
        randPoints.clear();

        RandomPointGenerator::apply(zb_it, x0, L_cov, Ntot, walk_length,
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


////////////////////////////////////
///
/// ratio estimation

// template <typename NT>
// bool is_max_error(NT const& a, NT const& b, NT const& error)
// {
//     return ((b-a)/a<error/2.0) ? true : false;
// }

// template <typename NT>
// struct estimate_ratio_parameters
// {
// public:

//     estimate_ratio_parameters(unsigned int W_len, unsigned int N, NT ratio)
//         :   min_val(std::numeric_limits<NT>::lowest())
//         ,   max_val(std::numeric_limits<NT>::max())
//         ,   max_iterations_estimation(10000000)
//         ,   min_index(W_len-1)
//         ,   max_index(W_len-1)
//         ,   W(W_len)
//         ,   index(0)
//         ,   tot_count(N)
//         ,   count_in(N * ratio)
//         ,   iter(0)
//         ,   last_W(std::vector<NT>(W_len))
//         ,   minmaxIt(last_W.begin())
//     {}

//     NT min_val;
//     NT max_val;
//     const unsigned int max_iterations_estimation;
//     unsigned int min_index;
//     unsigned int max_index;
//     unsigned int W;
//     unsigned int index;
//     size_t tot_count;
//     size_t count_in;
//     unsigned int iter;
//     std::vector<NT> last_W;
//     typename std::vector<NT>::iterator minmaxIt;
// };

// template <typename PolyEllipsoid, typename Point, typename NT>
// bool estimate_ratio_generic(PolyEllipsoid const& Pe2, Point const& p, NT const& error,
//        estimate_ratio_parameters<NT> &ratio_parameters)
// {
//     if (ratio_parameters.iter++ <= ratio_parameters.max_iterations_estimation)
//     {
//         if (Pe2.is_in(p) == -1) ratio_parameters.count_in = ratio_parameters.count_in + 1.0;

//         ratio_parameters.tot_count = ratio_parameters.tot_count + 1.0;
//         NT val = NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
//         ratio_parameters.last_W[ratio_parameters.index] = val;

//         if (val <= ratio_parameters.min_val)
//         {
//             ratio_parameters.min_val = val;
//             ratio_parameters.min_index = ratio_parameters.index;
//         } else if (ratio_parameters.min_index == ratio_parameters.index)
//         {
//             ratio_parameters.minmaxIt = std::min_element(ratio_parameters.last_W.begin(), ratio_parameters.last_W.end());
//             ratio_parameters.min_val = (*ratio_parameters.minmaxIt);
//             ratio_parameters.min_index = std::distance(ratio_parameters.last_W.begin(), ratio_parameters.minmaxIt);
//         }

//         if (val >= ratio_parameters.max_val)
//         {
//             ratio_parameters.max_val = val;
//             ratio_parameters.max_index = ratio_parameters.index;
//         } else if (ratio_parameters.max_index == ratio_parameters.index)
//         {
//             ratio_parameters.minmaxIt = std::max_element(ratio_parameters.last_W.begin(), ratio_parameters.last_W.end());
//             ratio_parameters.max_val = (*ratio_parameters.minmaxIt);
//             ratio_parameters.max_index = std::distance(ratio_parameters.last_W.begin(), ratio_parameters.minmaxIt);
//         }

//         if ( (ratio_parameters.max_val - ratio_parameters.min_val) / ratio_parameters.max_val <= error/2.0 )
//         {
//             return true;
//         }

//         ratio_parameters.index = ratio_parameters.index % ratio_parameters.W + 1;
//         if (ratio_parameters.index == ratio_parameters.W) ratio_parameters.index = 0;

//         return false;
//     }
//     return true;
// }


template
<
        typename WalkType,
        typename Point,
        typename PolyEllipsoid1,
        typename PolyEllipsoid2,
        typename MT,
        typename NT,
        typename RNG
>
NT estimate_ratio_ellipsoid(PolyEllipsoid1 const& Pe1,
                            PolyEllipsoid2 const& Pe2,
                            MT const& L_cov,
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
    WalkType walk(Pe1, p, L_cov, rng);

    do
    {
        walk.template apply(Pe1, p, L_cov, walk_length, rng);
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
    unsigned int n = E.dimension();
    Point p(n);

    typename Ellipsoid::VT eigenvals_inv_sqrt = E.eigenvals_inv_sqrt();
    typename Ellipsoid::MT eigenvecs = E.eigenvecs();

    do
    {
        p = GetPointInDellipsoid<Point>::apply(n, eigenvals_inv_sqrt, eigenvecs, rng);
    } while(!estimate_ratio_generic(Pe2, p, error, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}

//--------------------------------------------------------------------------

// template <typename NT>
// struct estimate_ratio_interval_parameters
// {
// public:
//     estimate_ratio_interval_parameters(unsigned int W_len,
//                                        unsigned int N,
//                                        NT ratio)
//         :    mean(0)
//         ,    sum_sq(0)
//         ,    sum(0)
//         ,    s(0)
//         ,    max_iterations_estimation(10000000)
//         ,    W(W_len)
//         ,    index(0)
//         ,    tot_count(N)
//         ,    count_in(N * ratio)
//         ,    iter(0)
//         ,    last_W(std::vector<NT>(W_len))
//     {}

//     NT mean;
//     NT sum_sq;
//     NT sum;
//     NT s;
//     const unsigned int max_iterations_estimation;
//     unsigned int W;
//     unsigned int index;
//     size_t tot_count;
//     size_t count_in;
//     unsigned int iter;
//     std::vector<NT> last_W;
// };

// template <typename PolyEllipsoid, typename Point, typename NT>
// void full_sliding_window(PolyEllipsoid const& Pe2,
//                          Point const& p,
//                          estimate_ratio_interval_parameters<NT>& ratio_parameters)
// {
//     if (Pe2.is_in(p) == -1) ratio_parameters.count_in = ratio_parameters.count_in + 1.0;

//     ratio_parameters.tot_count = ratio_parameters.tot_count + 1.0;
//     NT val = NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
//     ratio_parameters.sum += val;
//     ratio_parameters.sum_sq += val * val;
//     ratio_parameters.last_W[ratio_parameters.index] = val;
//     ratio_parameters.index = ratio_parameters.index % ratio_parameters.W + 1;
//     if (ratio_parameters.index == ratio_parameters.W) ratio_parameters.index = 0;
// }

// template <typename PolyEllipsoid, typename Point, typename NT>
// bool estimate_ratio_interval_generic(PolyEllipsoid const& Pe2,
//                                      Point const& p,
//                                      NT const& error,
//                                      NT const& zp,
//                                      estimate_ratio_interval_parameters
//                                      <NT>& ratio_parameters)
// {
//     if (ratio_parameters.iter++ <= ratio_parameters.max_iterations_estimation)
//     {
//         if (Pe2.is_in(p) == -1) ratio_parameters.count_in = ratio_parameters.count_in + 1.0;

//         ratio_parameters.tot_count = ratio_parameters.tot_count + 1.0;
//         NT val = NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);

//         ratio_parameters.mean = (ratio_parameters.mean
//                              - ratio_parameters.last_W[ratio_parameters.index] /
//                 NT(ratio_parameters.W)) + val / NT(ratio_parameters.W);

//         ratio_parameters.sum_sq = (ratio_parameters.sum_sq -
//                 ratio_parameters.last_W[ratio_parameters.index]
//                 * ratio_parameters.last_W[ratio_parameters.index])
//                 + val * val;

//         ratio_parameters.sum = (ratio_parameters.sum
//                                 - ratio_parameters.last_W[ratio_parameters.index])
//                                + val;

//         ratio_parameters.s = std::sqrt((ratio_parameters.sum_sq + NT(ratio_parameters.W) *
//                 ratio_parameters.mean * ratio_parameters.mean - NT(2)
//                                       * ratio_parameters.mean
//                                       * ratio_parameters.sum) /
//                                        NT(ratio_parameters.W));

//         ratio_parameters.last_W[ratio_parameters.index] = val;

//         ratio_parameters.index = ratio_parameters.index % ratio_parameters.W + 1;
//         if (ratio_parameters.index == ratio_parameters.W)
//         {
//             ratio_parameters.index = 0;
//         }

//         if (is_max_error(val - zp * ratio_parameters.s,
//                          val + zp * ratio_parameters.s,
//                          error))
//         {
//             return true;
//         }
//         return false;
//     }
//     return true;
// }

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

    // for sampling from ellipsoid
    typename Ellipsoid::VT eigenvals_inv_sqrt = E.eigenvals_inv_sqrt();
    typename Ellipsoid::MT eigenvecs = E.eigenvecs();

    unsigned int n = Pe2.dimension();
    Point p(n);

    for (int i = 0; i < ratio_parameters.W; ++i)
    {
        p = GetPointInDellipsoid<Point>::apply(n, eigenvals_inv_sqrt, eigenvecs, rng);
        full_sliding_window(Pe2, p, ratio_parameters);
    }

    ratio_parameters.mean = ratio_parameters.sum / NT(ratio_parameters.W);

    do {
        p = GetPointInDellipsoid<Point>::apply(n, eigenvals_inv_sqrt, eigenvecs, rng);
    } while (!estimate_ratio_interval_generic(Pe2, p, error, zp, ratio_parameters));

    return NT(ratio_parameters.count_in) / NT(ratio_parameters.tot_count);
}


template
<
        typename WalkType,
        typename Point,
        typename PolyEllipsoid1,
        typename PolyEllipsoid2,
        typename MT,
        typename NT,
        typename RNG
>
NT estimate_ratio_interval_ellipsoid(PolyEllipsoid1 const& Pe1,
                                     PolyEllipsoid2 const& Pe2,
                                     MT const& L_cov,
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
    WalkType walk(Pe1, p, L_cov, rng);

    for (int i = 0; i < ratio_parameters.W; ++i)
    {
        walk.template apply(Pe1, p, L_cov, walk_length, rng);
        full_sliding_window(Pe2, p, ratio_parameters);
    }

    ratio_parameters.mean = ratio_parameters.sum / NT(ratio_parameters.W);

    do {
        walk.template apply(Pe1, p, L_cov, walk_length, rng);
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
                                       unsigned int const& win_len = 300)
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
    VT x0 = P.inner_point();
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

    // ! std::pair<Point, NT> InnerBall = P.ComputeInnerBall();
    // ! if (InnerBall.second < 0.0) return std::pair<NT, NT> (-1.0, 0.0);

    // ! NT radius = InnerBall.second;
    // ! Point c = InnerBall.first;

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

    // ! NT vol = (NT(n)/NT(2) * std::log(M_PI)) + NT(n)*std::log((*(BallSet.end() - 1)).radius()) - log_gamma_function(NT(n) / NT(2) + 1);
    NT vol = (*(EllipsoidSet.end() - 1)).log_volume();

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

    auto ellipsoiditer = EllipsoidSet.begin();
    auto ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1)
    {
        vol += (!parameters.window2) ?
               std::log(NT(1) / estimate_ratio_interval_ellipsoid
                    <WalkType, Point>(P,
                                      *ellipsoiditer,
                                      L_cov,
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
                                      L_cov,
                                      *ratioiter,
                                      er1,
                                      parameters.win_len,
                                      N_times_nu,
                                      walk_length,
                                      rng));
    }

    for ( ; ellipsoiditer < EllipsoidSet.end() - 1; ++ellipsoiditer, ++ratioiter)
    {
        PolyEllipsoid Pe(P, *ellipsoiditer);
        vol += (!parameters.window2) ?
                    std::log(NT(1) / estimate_ratio_interval_ellipsoid
                                <WalkType, Point>(Pe,
                                                  *(ellipsoiditer + 1),
                                                  L_cov,
                                                  *(ratioiter + 1),
                                                  er1, parameters.win_len,
                                                  N_times_nu,
                                                  prob, walk_length,
                                                  rng))
                  : std::log(NT(1) / estimate_ratio_ellipsoid
                                <WalkType, Point>(Pe,
                                                  *ellipsoiditer,
                                                  L_cov,
                                                  *ratioiter,
                                                  er1,
                                                  parameters.win_len,
                                                  N_times_nu,
                                                  walk_length,
                                                  rng));
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
                                               unsigned int const& walk_length = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_cooling_ellipsoids<WalkTypePolicy>(Pin, rng, error, walk_length);
}


#endif // VOLUME_COOLING_ELLIPSOIDS_HPP
