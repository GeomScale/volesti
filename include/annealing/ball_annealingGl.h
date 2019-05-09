// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALINGGL_H
#define BALL_ANNEALINGGL_H

#include <boost/math/distributions/students_t.hpp>


template <class Point, class ConvexBody, class PointList, typename NT>
bool check_converg001(ConvexBody &P, PointList &randPoints, NT lb, NT ub, bool &too_few, NT &ratio,
                      int nu, bool precheck, bool lastball) {

    std::vector<NT> ratios;
    NT countsIn = 0.0;
    int m = randPoints.size()/nu;
    std::pair<NT,NT> mv;
    NT T, rs;

    NT alpha_check = 0.01;
    bool print = true;

    int i = 1, count_sets=0;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){

        if (P.is_in(*pit)==-1) countsIn += 1.0;
        if (i % m == 0) {
            ratios.push_back(countsIn/m);
            countsIn = 0.0;
            if (ratios.size()>1 && precheck) {
                boost::math::students_t dist(ratios.size() - 1);
                mv = getMeanVariance(ratios);
                ratio = mv.first;
                rs = std::sqrt(mv.second);
                T = rs * (boost::math::quantile(boost::math::complement(dist, alpha_check / 2.0))
                          / std::sqrt(NT(ratios.size())));
                if (ratio + T < lb) {
                    too_few = true;
                    return false;
                } else if (ratio - T > ub) return false;
            }
        }
    }

    NT alpha = 0.25;
    if(precheck) alpha = 0.1;
    mv = getMeanVariance(ratios);
    ratio = mv.first;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs*(boost::math::quantile(boost::math::complement(dist, alpha))
              / std::sqrt(NT(nu)));
    if (ratio > lb + T) {
        if (lastball) return true;
        if ((precheck && ratio < ub - T) || (!precheck && ratio < ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;

}


template <class Point, class ball, class PointList, typename NT>
void get_next_zonoball(std::vector<ball> &BallSet,
                       PointList &randPoints, NT rad_min, std::vector<NT> &ratios, NT lb, NT ub, int nu){

    int n = (*randPoints.begin()).dimension();
    bool too_few;
    NT radmax = 0.0, rad, pnorm, ratio;

    for (typename PointList::iterator rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit) {
        pnorm = (*rpit).squared_length();
        if (pnorm > radmax) radmax = pnorm;
    }
    ball Biter;
    radmax=std::sqrt(radmax);

    while (true) {
        rad = 0.5 * (rad_min + radmax);
        Biter = ball(Point(n), rad * rad);
        too_few = false;

        if (check_converg001<Point>(Biter, randPoints, lb, ub, too_few, ratio, nu, false, false)) {
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return;
        }

        if (too_few) {
            rad_min = rad;
        } else {
            radmax = rad;
        }
    }

}

template <class RNGType,class ball, class Polytope, typename NT>
void get_first_ball(Polytope &P, ball &B0, NT &ratio, NT radius, NT lb, NT ub, NT rmax){

    typedef typename Polytope::PolytopePoint Point;
    int n = P.dimension();
    bool bisection_int = false, pass = false, too_few = false;
    bool print = true;
    std::list<Point> randPoints;
    Point p(n);

    if(rmax>0.0) {
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));
        }
        pass = check_converg001<Point>(P, randPoints, lb, ub, too_few, ratio, 10, true, false);
        if (pass || !too_few) {
            B0 = ball(Point(n), rmax*rmax);
            return;
        }
        bisection_int = true;
    } else {
        rmax = 2 * std::sqrt(NT(n)) * radius;
    }
    NT rad1 = radius;

    while(!bisection_int) {

        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));

        if(check_converg001<Point>(P, randPoints, lb, ub, too_few, ratio, 10, true, false)) {
            B0 = ball(Point(n), rmax*rmax);
            return;
        }

        if (too_few) break;
        rad1 = rmax;
        rmax = rmax + 2*std::sqrt(NT(n))*radius;
    }

    NT rad_med;

    while(true) {

        rad_med = 0.5*(rad1+rmax);
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rad_med));

        if(check_converg001<Point>(P, randPoints, lb, ub, too_few, ratio, 10, true, false)) {
            B0 = ball(Point(n), rad_med*rad_med);
            return;
        }

        if (too_few) {
            rmax = rad_med;
        } else {
            rad1 = rad_med;
        }

    }

}

template <class PolyBall, class RNGType,class ball, class Polytope, class Parameters, typename NT>
void get_sequence_of_zonoballs(Polytope &P, std::vector<ball> &BallSet, ball &B0, NT &ratio0,
                               std::vector<NT> &ratios, int Ntot, int nu,
                               NT lb, NT ub, NT radius, Parameters &var,
                               NT rmax = 0.0) {


    typedef typename Polytope::PolytopePoint Point;
    typedef typename Polytope::MT MT;
    bool print = var.verbose, pass, fail;
    int n = P.dimension();
    NT ratio;
    std::list<Point> randPoints;
    Point q(n);
    PolyBall zb_it;
    get_first_ball<RNGType>(P, B0, ratio, radius, lb, ub, rmax);
    ratio0 = ratio;
    rand_point_generator(P, q, Ntot, var.walk_steps, randPoints, var);

    if (check_converg001<Point>(B0, randPoints, lb, ub, fail, ratio, nu, false, true)) {
        ratios.push_back(ratio);
        if(print) std::cout<<"one ball and ratio = "<<ratio<<std::endl;
        return;
    }
    if(print) std::cout<<"not the last ball, ratio = "<<ratio<<std::endl;
    get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, nu);
    if(print) std::cout<<"number of balls = "<<BallSet.size()<<std::endl;

    while (true) {
        zb_it = PolyBall(P, BallSet[BallSet.size()-1]);
        q=Point(n);
        randPoints.clear();
        rand_point_generator(zb_it, q, Ntot, var.walk_steps, randPoints, var);

        if (check_converg001<Point>(B0, randPoints, lb, ub, fail, ratio, nu, false, true)) {
            ratios.push_back(ratio);
            return;
        }
        get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, nu);
        if(print) std::cout<<"number of balls = "<<BallSet.size()<<std::endl;
    }

}


#endif
