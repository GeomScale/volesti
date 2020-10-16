// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_COOLING_HPOLY_HPP
#define VOLUME_COOLING_HPOLY_HPP

#include "volume/volume_cooling_gaussians.hpp"
#include "sampling/random_point_generators.hpp"


template
<
        typename RandomPointGenerator,
        typename Zonotope,
        typename HPolytope,
        typename NT,
        typename RNG,
        typename VT
>
bool get_first_poly(Zonotope &P,
                    HPolytope &HP,
                    NT &ratio,
                    cooling_ball_parameters<NT> const& parameters,
                    RNG& rng, VT &b_max)
{

    typedef typename Zonotope::PointType Point;
    typedef typename Zonotope::MT MT;
    //typedef typename Zonotope::VT VT;

    PushBackWalkPolicy push_back_policy;

    const unsigned max_iterarions = 20;
    NT tolerance = 0.00000000001;

    MT G = P.get_mat().transpose(), A = HP.get_mat();
    b_max = (A*G).cwiseAbs().rowwise().sum();
    VT b_min = HP.get_vec();
    HPolytope HPiter=HP;

    int n = P.dimension(), m = b_max.size(), N = 1200, iter = 1, count = 0;
    Point q(n);
    bool too_few, print = false;
    std::list<Point> randPoints;

    NT l=0.0, u=1.0, med;
    VT  b_med(m);

    while(iter <= max_iterarions) {

        q=Point(n);
        med = (u + l) * 0.5;
        b_med = b_min + (b_max-b_min)*med;
        HPiter.set_vec(b_med);

        randPoints.clear();
        RandomPointGenerator::apply(HPiter, q, 1200, 10+10*n,
                                    randPoints, push_back_policy, rng);
        //rand_point_generator(HPiter, q, 1200, 10+10*n, randPoints, var);
        too_few = false;

        if(check_convergence<Point>(P, randPoints,
                                    too_few, ratio, parameters.nu,
                                 true, false, parameters)){
        //if(check_convergence<Point>(P, randPoints, lb, up_lim, too_few, ratio, 10, 0.2, true, false)) {
            HP.set_vec(b_med);
            return true;
        }

        if (too_few) {
            u = med;
        } else {
            l = med;
        }
        if(med>0.9) {
            HP.set_vec(b_med);
            return true;
        }
        if(u-l < tolerance) {
            u=1.0;
            l=0.0;
            iter++;
        }
    }
    return false;
}


template <typename Zonotope, typename HPolytope, typename VT, typename PointList, typename NT>
bool get_next_zonoball(std::vector<HPolytope> &HPolySet,
                       HPolytope &HP2,
                       const VT &b_max,
                       const VT &b_min,
                       PointList &randPoints,
                       std::vector<NT> &ratios,
                       cooling_ball_parameters<NT> const& parameters)
{

    typedef typename Zonotope::PointType Point;

    const unsigned max_iterarions = 20;
    NT tolerance = 0.00000000001;

    int n = HP2.dimension(), iter = 1;
    bool too_few;
    VT b_med(b_max.size());
    NT ratio, med, u = 1.0, l = 0.0;

    while (iter <= max_iterarions) {
        med = (u + l) * 0.5;
        b_med = b_min + (b_max-b_min) * med;
        HP2.set_vec(b_med);
        too_few = false;

        if(check_convergence<Point>(HP2, randPoints,
                                    too_few, ratio, parameters.nu,
                                    false, false, parameters)){
            HPolySet.push_back(HP2);
            ratios.push_back(ratio);
            return true;
        }
        if(too_few) {
            l = med;
        } else {
            u = med;
        }
        if(u-l < tolerance) {
            u=1.0;
            l=0.0;
            iter++;
        }
    }
    return false;
}

template <
        typename RandomPointGenerator,
        typename ZonoHP,
        typename Zonotope,
        typename HPolytope,
        typename VT,
        typename NT,
        typename RNG>
bool get_sequence_of_zonopolys(Zonotope &Z,
                               const HPolytope &HP,
                               std::vector<HPolytope> &HPolySet,
                               std::vector<NT> &ratios,
                               const VT &b_max,
                               unsigned int const& N_times_nu,
                               unsigned int const& walk_length,
                               cooling_ball_parameters<NT> const& parameters,
                               RNG& rng)
{

    bool too_few=false;
    typedef typename Zonotope::PointType Point;
    typedef typename Zonotope::MT MT;

    PushBackWalkPolicy push_back_policy;

    int n = Z.dimension();
    MT G = Z.get_mat().transpose();
    MT AG = HP.get_mat()*G;
    NT ratio;
    std::list<Point> randPoints;
    Point q(n);

    RandomPointGenerator::apply(Z, q, N_times_nu, walk_length,
                                randPoints, push_back_policy, rng);
    HPolytope HP2 = HP;
    if (check_convergence<Point>(HP, randPoints,
                                 too_few, ratio, parameters.nu,
                                 false, true, parameters)) {
        ratios.push_back(ratio);
        return true;
    }
    if ( !get_next_zonoball<Zonotope>(HPolySet, HP2, b_max, HP.get_vec(),
                                      randPoints, ratios, parameters))
    {
        return false;
    }

    ZonoHP ZHP2;
    VT Zs_min = HP.get_vec();

    while (true) {

        ZHP2 = ZonoHP(Z,HP2);
        q=Point(n);
        randPoints.clear();
        RandomPointGenerator::apply(ZHP2, q, N_times_nu, walk_length,
                                    randPoints, push_back_policy, rng);

        if (check_convergence<Point>(HP, randPoints,
                                     too_few, ratio, parameters.nu,
                                     false, true, parameters))
        {
            ratios.push_back(ratio);
            return true;
        }
        if ( !get_next_zonoball<Zonotope>(HPolySet, HP2, HP2.get_vec(),
                                          Zs_min, randPoints, ratios, parameters) )
        {
            return false;
        }
    }
}


template
        <
                typename Zonotope,
                typename HPolytope
        >
void compute_hpoly_for_mmc(Zonotope &P, HPolytope &HP) {

    typedef typename Zonotope::PointType Point;
    typedef typename Zonotope::NT NT;
    typedef typename Zonotope::VT VT;
    typedef typename Zonotope::MT MT;

    MT V = P.get_mat();
    MT G = V.transpose();
    int m = G.cols();
    std::list<Point> randPoints;

    MT XX(m, 2*m);
    XX << MT::Identity(m,m), -MT::Identity(m,m);
    MT AA = XX.transpose(); VT b = VT::Ones(2*m);
    MT T = P.get_T();
    MT Tt = T.transpose();
    MT A2 = AA * Tt, B = G * Tt;
    MT A3 = A2 * B.inverse();

    NT row_norm;
    for (int i = 0; i < A3.rows(); ++i) {
        row_norm = A3.row(i).norm();
        A3.row(i) = A3.row(i) / row_norm;
        b(i) = b(i) / row_norm;
    }

    HP.init(P.dimension(),A3,b);
}


template
        <
                typename WalkTypePolicy,
                typename HPolytope,
                typename Zonotope,
                typename RandomNumberGenerator
        >
double volume_cooling_hpoly (Zonotope const& Pin,
                         RandomNumberGenerator &rng,
                         double const& error = 0.1,
                         unsigned int const& walk_length = 1,
                             unsigned int const& win_len = 250)
{

    typedef typename Zonotope::PointType Point;
    typedef typename Point::FT NT;
    typedef ZonoIntersectHPoly <Zonotope, HPolytope> ZonoHP;
    typedef typename Zonotope::VT VT;
    typedef typename Zonotope::MT MT;
    typedef std::list <Point> PointList;

    typedef typename WalkTypePolicy::template Walk
            <
                    Zonotope,
                    RandomNumberGenerator
            > WalkType;
    typedef RandomPointGenerator<WalkType> ZonoRandomPointGenerator;

    typedef typename CDHRWalk::template Walk
            <
                    HPolytope,
                    RandomNumberGenerator
            > CdhrWalk;
    typedef RandomPointGenerator<CdhrWalk> CdhrRandomPointGenerator;

    auto P(Pin);
    //RandomNumberGenerator rng(P.dimension());
    //RandomNumberGenerator rng_diam(P.num_of_generators());
    cooling_ball_parameters<NT> parameters(win_len);

    int n = P.dimension();
    NT prob = parameters.p, ratio;
    int N_times_nu = parameters.N * parameters.nu;

    HPolytope HP;
    compute_hpoly_for_mmc(P, HP);
    VT b_max(2*P.num_of_generators());
    if ( !get_first_poly<CdhrRandomPointGenerator>(P, HP, ratio, parameters, rng, b_max) )
    {
        return -1.0;
    }

    std::vector<HPolytope > HPolySet;
    std::vector<NT> ratios;

    ZonoHP zb1, zb2;
    std::vector<NT> diams_inter;

    if ( !get_sequence_of_zonopolys<ZonoRandomPointGenerator, ZonoHP>
                       (P, HP, HPolySet, ratios,
                        b_max, N_times_nu, walk_length, parameters, rng) )
    {
        return -1.0;
    }

    int mm=HPolySet.size()+2;
    int mm2=mm+1;
    prob = std::pow(prob, 1.0/NT(mm2));
    NT er0 = error/(2.0*std::sqrt(NT(mm2)));
    NT er1 = (error*std::sqrt(2.0*NT(mm2)-1))/(std::sqrt(2.0*NT(mm2)));
    NT Her = error/(2.0*std::sqrt(NT(mm2)));


    HPolytope HP2(HP);
    std::pair<Point, NT> InnerBall = HP2.ComputeInnerBall();
    std::pair< std::pair<MT, VT>, NT > res = round_polytope<CDHRWalk, MT, VT>(HP2, InnerBall,
            10 + 10 * n, rng);
    //TODO: rounding to HP2
    //NT vol = res.second * volume_cooling_balls<BilliardWalk, RandomNumberGenerator>(HP2, Her, 1);
    NT vol = res.second * volume_cooling_gaussians<GaussianCDHRWalk>(HP2, rng, Her/2.0, 1);

    if (!parameters.window2) {
        vol *= estimate_ratio_interval<CdhrWalk, Point>(HP, P, ratio, er0, parameters.win_len, 1200, prob, 10+10*n, rng);
    } else {
        vol *= estimate_ratio<CdhrWalk, Point>(HP, P, ratio, er0, parameters.win_len, 1200, 10 + 10*n, rng);
    }

    HPolytope b1, b2;
    if (HPolySet.size()==0) {
        if (ratios[0]!=1) {
            if(!parameters.window2) {
                vol = vol / estimate_ratio_interval<WalkType, Point>(P, HP, ratios[0], er1, parameters.win_len, N_times_nu,
                                                               prob, walk_length, rng);
            } else {
                vol = vol / estimate_ratio<WalkType, Point>(P, HP, ratios[0], er1, parameters.win_len, N_times_nu,
                                                               walk_length, rng);
            }
        }
    } else {
        er1 = er1 / std::sqrt(NT(mm)-1.0);
        b1 = HPolySet[0];
        if(!parameters.window2) {
            vol = vol / estimate_ratio_interval<WalkType, Point>(P, b1, ratios[0], er1, parameters.win_len, N_times_nu,
                                                               prob, walk_length, rng);
        } else {
            vol = vol / estimate_ratio<WalkType, Point>(P, b1, ratios[0], er1, parameters.win_len, N_times_nu,
                                                               walk_length, rng);
        }

        for (int i = 0; i < HPolySet.size()-1; ++i) {
            zb1 = ZonoHP(P,HPolySet[i]);
            b2 = HPolySet[i+1];
            if(!parameters.window2) {
                vol = vol / estimate_ratio_interval<WalkType, Point>(zb1, b2, ratios[i], er1, parameters.win_len,
                                                               N_times_nu, prob, walk_length, rng);
            } else {
                vol = vol / estimate_ratio<WalkType, Point>(zb1, b2, ratios[i], er1, parameters.win_len, N_times_nu,
                                                               walk_length, rng);
            }
        }

        zb1 = ZonoHP(P,HPolySet[HPolySet.size()-1]);
        if (!parameters.window2) {
            vol = vol / estimate_ratio_interval<WalkType, Point>(zb1, HP, ratios[ratios.size() - 1], er1,
                                                               parameters.win_len, N_times_nu, prob, walk_length, rng);
        } else {
            vol = vol / estimate_ratio<WalkType, Point>(zb1, HP, ratios[ratios.size() - 1], er1, parameters.win_len,
                                                               N_times_nu, walk_length, rng);
        }
    }

    P.free_them_all();

    return vol;

}


template
<
        typename WalkTypePolicy,
        typename RandomNumberGenerator,
        typename HPolytope,
        typename Polytope
>
double volume_cooling_hpoly(Polytope const& Pin,
                                 double const& error = 0.1,
                                 unsigned int const& walk_length = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_cooling_hpoly<WalkTypePolicy, HPolytope>(Pin, rng, error, walk_length);
}

#endif // VOLUME_COOLING_HPOLY_HPP
