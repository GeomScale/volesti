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

template <typename Point, typename ConvexBody, typename PointList, typename NT>
bool check_convergence2(ConvexBody &P,
                        PointList &randPoints,
                        const NT &lb,
                        const NT &ub,
                        bool &too_few,
                        NT &ratio,
                        const int &nu,
                        NT alpha,
                        const bool &precheck,
                        const bool &lastball)
{

    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int m = randPoints.size()/nu, i = 1;
    NT T, rs, alpha_check = 0.01;
    size_t countsIn = 0;

    for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++)
    {

        if (P.is_in(*pit)==-1) countsIn++;
        if (i % m == 0) {
            ratios.push_back(NT(countsIn)/m);
            countsIn = 0;
            if (ratios.size()>1 && precheck) {
                boost::math::students_t dist(ratios.size() - 1);
                mv = getMeanVariance(ratios);
                ratio = mv.first;
                rs = std::sqrt(mv.second);
                T = rs * (boost::math::quantile
                            (boost::math::complement(dist, alpha_check / 2.0))
                          / std::sqrt(NT(ratios.size())));
                if (ratio + T < lb) {
                    too_few = true;
                    return false;
                } else if (ratio - T > ub) return false;
            }
        }
    }

    //NT alpha = 0.25;
    if(precheck) alpha *= 0.5;
    mv = getMeanVariance(ratios);
    ratio = mv.first;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs*(boost::math::quantile(boost::math::complement(dist, alpha))
            / std::sqrt(NT(nu)));
    if (ratio > lb + T)
    {
        if (lastball) return true;
        if ((precheck && ratio < ub - T) || (!precheck && ratio < ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;
}


template <typename Polytope, typename ball, typename NT, typename RNG>
bool get_first_ball2(Polytope &P,
                     ball &B0,
                     NT &ratio,
                     NT rad1,
                     const NT &lb,
                     const NT &ub,
                     const NT &alpha,
                     NT &rmax,
                     RNG& rng)
{

    typedef typename Polytope::PointType Point;
    int n = P.dimension(), iter = 1;
    bool bisection_int = false, pass = false, too_few = false;
    std::list<Point> randPoints;
    Point p(n);

    if (rmax>0.0)
    {
        for (int i = 0; i < 1200; ++i)
        {
            //randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rmax, rng));
        }
        pass = check_convergence2<Point>(P, randPoints, lb, ub, too_few, ratio,
                                         10, alpha, true, false);
        if (pass || !too_few) {
            B0 = ball(Point(n), rmax*rmax);
            return true;
        }
        bisection_int = true;
    } else {
        rmax = 2 * std::sqrt(NT(n)) * rad1;
    }
    NT radius = rad1;

    while(!bisection_int)
    {

        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i)
        {
            //randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rmax, rng));
        }

        if (check_convergence2<Point>(P, randPoints, lb, ub, too_few, ratio, 10,
                                      alpha, true, false))
        {
            B0 = ball(Point(n), rmax*rmax);
            return true;
        }

        if (too_few) break;
        rad1 = rmax;
        rmax = rmax + 2*std::sqrt(NT(n))*radius;
    }

    NT rad_med, rad0=rad1, rad_m = rmax;

    while(iter <= MAX_ITER)
    {

        rad_med = 0.5*(rad1+rmax);
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i)
        {
            //randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rad_med));
            randPoints.push_back(GetPointInDsphere<Point>::apply(n, rad_med, rng));
        }

        if (check_convergence2<Point>(P, randPoints, lb, ub, too_few, ratio, 10,
                                      alpha, true, false))
        {
            B0 = ball(Point(n), rad_med*rad_med);
            return true;
        }

        if (too_few) {
            rmax = rad_med;
        } else {
            rad1 = rad_med;
        }

        if(rmax-rad1 < TOL) {
            rad1 = rad0;
            rmax = rad_m;
            iter++;
        }

    }
    return false;
}


template <typename Point, typename ball, typename PointList, typename NT>
bool get_next_zonotopeball(std::vector<ball> &BallSet,
                           PointList &randPoints,
                           NT rad_min,
                           std::vector<NT> &ratios,
                           const NT &lb,
                           const NT &ub,
                           NT &alpha,
                           const int &nu)
{

    int n = (*randPoints.begin()).dimension(), iter = 1;
    bool too_few;
    NT radmax = 0.0, rad, pnorm, ratio;

    for (auto rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit)
    {
        pnorm = (*rpit).squared_length();
        if (pnorm > radmax) radmax = pnorm;
    }
    ball Biter;
    radmax=std::sqrt(radmax);
    NT rad0 = rad_min, rad_m = radmax;

    while (iter <= MAX_ITER) {
        rad = 0.5 * (rad_min + radmax);
        Biter = ball(Point(n), rad * rad);
        too_few = false;

        if (check_convergence2<Point>(Biter, randPoints, lb, ub, too_few, ratio,
                                      nu, alpha, false, false))
        {
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return true;
        }

        if (too_few) {
            rad_min = rad;
        } else {
            radmax = rad;
        }

        if(radmax-rad_min < TOL) {
            rad_min = rad0;
            radmax = rad_m;
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
bool get_sequence_of_polytopeballs(Polytope &P,
                                   std::vector<ball> &BallSet,
                                   std::vector<NT> &ratios,
                                   const int &Ntot,
                                   const int &nu,
                                   const NT &lb,
                                   const NT &ub,
                                   NT radius,
                                   NT &alpha,
                                   unsigned int const& walk_length,
                                   NT& diameter,
                                   NT &rmax,
                                   RNG& rng)
{

    typedef typename Polytope::PointType Point;
    bool fail;
    int n = P.dimension();
    NT ratio, ratio0;
    std::list<Point> randPoints;
    ball B0;
    Point q(n);
    PolyBall zb_it;

    if ( !get_first_ball2(P, B0, ratio, radius, lb, ub, alpha, rmax, rng) )
    {
        return false;
    }

    ratio0 = ratio;
    //rand_point_generator(P, q, Ntot, var.walk_steps, randPoints, var);

    PushBackWalkPolicy push_back_policy;
    RandomPointGenerator::apply(P, q, Ntot, walk_length,
                                randPoints, push_back_policy, rng);

    if (check_convergence2<Point>(B0, randPoints, lb, ub, fail, ratio, nu,
                                  alpha, false, true))
    {
        ratios.push_back(ratio);
        BallSet.push_back(B0);
        ratios.push_back(ratio0);
        return true;
    }
    if ( !get_next_zonotopeball<Point>(BallSet, randPoints, B0.radius(), ratios,
                                       lb, ub, alpha, nu) )
    {
        return false;
    }

    while (true)
    {
        zb_it = PolyBall(P, BallSet[BallSet.size()-1]);
        q=Point(n);
        randPoints.clear();
        zb_it.comp_diam(diameter);

        //rand_point_generator(zb_it, q, Ntot, var.walk_steps, randPoints,var);
        RandomPointGenerator::apply(zb_it, q, Ntot, walk_length,
                                    randPoints, push_back_policy, rng);
        if (check_convergence2<Point>(B0, randPoints, lb, ub, fail, ratio, nu,
                                      alpha, false, true))
        {
            ratios.push_back(ratio);
            BallSet.push_back(B0);
            ratios.push_back(ratio0);
            return true;
        }
        if ( !get_next_zonotopeball<Point>(BallSet, randPoints, B0.radius(),
                                           ratios, lb, ub, alpha, nu) )
        {
            return false;
        }
    }
}


////////////////////////////////////
///
/// ratio estimation

template <typename NT>
bool check_max_error2(const NT &a, const NT &b, const NT &error)
{

    if((b-a)/a<error/2.0) {
        return true;
    }
    return false;

}


template
<
    typename Point,
    typename PolyBall1,
    typename PolyBall2,
    typename NT,
    typename WalkType,
    typename RNG
>
NT estimate_ratio(PolyBall1 &Pb1,
                  PolyBall2 &Pb2,
                  const NT &ratio,
                  const NT &error,
                  const int &W,
                  const int &Ntot,
                  const unsigned int& walk_length,
                  WalkType& walk,
                  RNG& rng,
                  bool isball = false,
                  NT radius = 0.0)
{

    int n = Pb1.dimension();
    int min_index = W-1;
    int max_index = W-1;
    int index = 0;
    int iter = 1;
    NT min_val = std::numeric_limits<NT>::lowest();
    NT max_val = std::numeric_limits<NT>::max();
    NT val;
    size_t totCount = Ntot, countIn = Ntot * ratio;
    std::vector<NT> last_W(W);

    typename std::vector<NT>::iterator minmaxIt;
    Point p(n);

    //if (!var.ball_walk && !isball)
    //{
    //    uniform_first_point(Pb1,p,p_prev,coord_prev,var.walk_steps,
    //                        lamdas,Av,lambda,var);
    //}

    while(iter <= MAX_ITER_ESTI)
    {
        iter++;

        if (isball)
        {
            //p = get_point_in_Dsphere<RNGType, Point>(n, radius);
            p = GetPointInDsphere<Point>::apply(n, radius, rng);
        } else {
            //uniform_next_point(Pb1, p, p_prev, coord_prev, var.walk_steps,
            //                   lamdas, Av, lambda, var);
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
        } else if(min_index==index)
        {
            minmaxIt = std::min_element(last_W.begin(), last_W.end());
            min_val = *minmaxIt;
            min_index = std::distance(last_W.begin(), minmaxIt);
        }

        if (val>=max_val)
        {
            max_val = val;
            max_index = index;
        }else if (max_index==index)
        {
            minmaxIt = std::max_element(last_W.begin(), last_W.end());
            max_val = *minmaxIt;
            max_index = std::distance(last_W.begin(), minmaxIt);
        }

        if ( (max_val-min_val)/max_val<=error/2.0 )
        {
            return val;
        }

        index = index%W+1;
        if (index==W) index=0;
    }
    return val;
}


template
<
    typename Point,
    typename PolyBall1,
    typename PolyBall2,
    typename NT,
    typename WalkType,
    typename RNG
>
NT estimate_ratio_interval(PolyBall1 &Pb1,
                           PolyBall2 &Pb2,
                           const NT &ratio,
                           const NT &error,
                           const int &W,
                           const int &Ntot,
                           const NT &prob,
                           const unsigned int& walk_length,
                           WalkType& walk,
                           RNG& rng,
                           bool isball = false,
                           NT radius = 0.0)
{

    int n = Pb1.dimension();
    int index = 0;
    int iter = 1;
    std::vector<NT> last_W(W);
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    VT lamdas, Av;
    Av.setZero(Pb1.num_of_hyperplanes());
    lamdas.setZero(Pb1.num_of_hyperplanes());
    NT val, sum_sq=0.0, sum=0.0, lambda;
    size_t totCount = Ntot, countIn = Ntot * ratio;
    //std::cout<<"countIn = "<<countIn<<", totCount = "<<totCount<<std::endl;

    Point p(n);
    Point p_prev=p;
    unsigned int coord_prev;
    //if (!var.ball_walk && !isball) uniform_first_point(Pb1, p, p_prev, coord_prev, 1,
    //                                                  lamdas, Av, lambda, var);

    for (int i = 0; i < W; ++i)
    {
        if (isball)
        {
            //p = get_point_in_Dsphere<RNGType, Point>(n, radius);
            p = GetPointInDsphere<Point>::apply(n, radius, rng);
        } else {
            //uniform_next_point(Pb1, p, p_prev, coord_prev, var.walk_steps,
            //                   lamdas, Av, lambda, var);
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

    while(iter <= MAX_ITER_ESTI) {
        iter++;

        if (isball) {
            //p = get_point_in_Dsphere<RNGType, Point>(n, radius);
            p = GetPointInDsphere<Point>::apply(n, radius, rng);
        } else {
            //uniform_next_point(Pb1, p, p_prev, coord_prev, var.walk_steps,
            //                   lamdas, Av, lambda, var);
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
        if (check_max_error2(val - zp * s, val + zp * s, error)) {
            //if (print) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
            return val;
        }

    }
    return val;

}


template
<
    typename WalkTypePolicy = BilliardWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, double>,
    typename Polytope,
    typename AParameters
>
double volume_cooling_balls(Polytope &P,
                            AParameters &var_ban,
                            double const& error = 1.0,
                            unsigned int const& walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> ball;
    typedef BallIntersectPolytope <Polytope, ball> PolyBall;
    typedef typename Polytope::VT VT;
    typedef std::list <Point> PointList;

    typedef typename WalkTypePolicy::template Walk
                                              <
                                                Polytope,
                                                RandomNumberGenerator
                                              > WalkType;
    typedef RandomPointGenerator<WalkType> RandomPointGenerator;

    RandomNumberGenerator rng(P.dimension());
    int n = P.dimension();
    int win_len = var_ban.win_len;
    int N = var_ban.N;
    int nu = var_ban.nu;
    bool window2 = var_ban.window2;
    NT lb = var_ban.lb;
    NT ub = var_ban.ub;
    NT prob = var_ban.p;
    NT rmax = var_ban.rmax;
    auto InnerBall = P.InnerBall();
    NT radius = InnerBall.second;
    NT e = error;
    NT alpha = var_ban.alpha;
    NT diameter = P.ComputeDiameter();

    std::vector <ball> BallSet;
    std::vector <NT> ratios;
    Point c = InnerBall.first;
    P.normalize();

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(c.getCoefficients());

    WalkType walk(P, c, rng);

    if ( !get_sequence_of_polytopeballs
          <
            RandomPointGenerator,
            PolyBall
         >(P, BallSet, ratios,
           N * nu, nu, lb, ub,
           radius, alpha, walk_length, diameter, rmax, rng) )
    {
        return -1.0;
    }
    //var.diameter = diam;

    NT vol = (std::pow(M_PI, n / 2.0)
              * (std::pow((*(BallSet.end() - 1)).radius(), n)))
            / (tgamma(n / 2.0 + 1));

    int mm = BallSet.size() + 1;
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = e / (2.0 * std::sqrt(NT(mm)));
    NT er1 = (e * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    vol *= (window2) ?
           estimate_ratio<Point>(*(BallSet.end() - 1),
                                      P,
                                      *(ratios.end() - 1),
                                      er0, win_len, 1200, walk_length,
                                      walk, rng,
                                      true,
                                      (*(BallSet.end() - 1)).radius())
        :
           estimate_ratio_interval<Point>(*(BallSet.end() - 1),
                                               P,
                                               *(ratios.end() - 1),
                                               er0, win_len, 1200, prob,
                                               walk_length, walk, rng,
                                               true,
                                               (*(BallSet.end() - 1)).radius());

    PolyBall Pb;
    typename std::vector<ball>::iterator balliter = BallSet.begin();
    typename std::vector<NT>::iterator ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1)
    {
        vol *= (!window2) ?
               1 / estimate_ratio_interval<Point>(P, *balliter, *ratioiter,
                                                       er1, win_len, N * nu, prob,
                                                       walk_length, walk, rng)
            : 1 / estimate_ratio<Point>(P, *balliter, *ratioiter,
                                             er1, win_len, N * nu, walk_length, walk, rng);
    }
    for ( ; balliter < BallSet.end() - 1; ++balliter, ++ratioiter)
    {
        Pb = PolyBall(P, *balliter);
        Pb.comp_diam(diameter);
        vol *= (!window2) ?
               1 / estimate_ratio_interval<Point>(Pb,
                                                       *(balliter + 1),
                                                       *(ratioiter + 1),
                                                       er1, win_len, N * nu,
                                                       prob, walk_length, walk, rng)
            : 1 / estimate_ratio<Point>(Pb, *balliter, *ratioiter, er1,
                                             win_len, N * nu, walk_length, walk, rng);
    }

    P.free_them_all();
    return vol;
}

#endif
