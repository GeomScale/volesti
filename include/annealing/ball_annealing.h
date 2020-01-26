// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALING_H
#define BALL_ANNEALING_H

#define MAX_ITER 10
#define TOL 0.00000000001


template <class Point, class ConvexBody, class PointList, typename NT>
bool check_convergence(ConvexBody &P, PointList &randPoints, const NT &lb, const NT &ub, bool &too_few, NT &ratio,
                      const int &nu, NT alpha, const bool &precheck, const bool &lastball) {

    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int m = randPoints.size()/nu;
    NT T, rs, alpha_check = 0.01, countsIn = 0.0;

    bool print = true;

    int i = 1;
    for(typename PointList::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){

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

    //NT alpha = 0.25;
    if(precheck) alpha *= 0.5;
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
bool get_next_zonoball(std::vector<ball> &BallSet, PointList &randPoints, NT rad_min, std::vector<NT> &ratios,
                       const NT &lb, const NT &ub, NT &alpha, const int &nu){

    int n = (*randPoints.begin()).dimension(), iter = 1;
    bool too_few;
    NT radmax = 0.0, rad, pnorm, ratio;

    for (typename PointList::iterator rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit) {
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

        if (check_convergence<Point>(Biter, randPoints, lb, ub, too_few, ratio, nu, alpha, false, false)) {
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


template <class RNGType, class Polytope, class ball, typename NT>
bool get_first_ball(Polytope &P, ball &B0, NT &ratio, NT rad1, const NT &lb, const NT &ub, const NT &alpha, NT &rmax){

    typedef typename Polytope::PolytopePoint Point;
    int n = P.dimension(), iter = 1;
    bool bisection_int = false, pass = false, too_few = false;
    std::list<Point> randPoints;
    Point p(n);

    if(rmax>0.0) {
        for (int i = 0; i < 1200; ++i) {
            randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));
        }
        pass = check_convergence<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false);
        if (pass || !too_few) {
            B0 = ball(Point(n), rmax*rmax);
            return true;
        }
        bisection_int = true;
    } else {
        rmax = 2 * std::sqrt(NT(n)) * rad1;
    }
    NT radius = rad1;

    while(!bisection_int) {

        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rmax));

        if(check_convergence<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false)) {
            B0 = ball(Point(n), rmax*rmax);
            return true;
        }

        if (too_few) break;
        rad1 = rmax;
        rmax = rmax + 2*std::sqrt(NT(n))*radius;
    }

    NT rad_med, rad0=rad1, rad_m = rmax;

    while(iter <= MAX_ITER) {

        rad_med = 0.5*(rad1+rmax);
        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i) randPoints.push_back(get_point_in_Dsphere<RNGType, Point>(n, rad_med));

        if(check_convergence<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false)) {
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

template <class PolyBall, class RNGType,class ball, class Polytope, class Parameters, typename NT>
bool get_sequence_of_polyballs(Polytope &P, std::vector<ball> &BallSet, std::vector<NT> &ratios, const int &Ntot, const int &nu,
                               const NT &lb, const NT &ub, NT radius, NT &alpha, Parameters &var, NT &rmax) {

    typedef typename Polytope::PolytopePoint Point;
    typedef typename Polytope::MT MT;
    bool fail;
    int n = P.dimension();
    NT ratio, ratio0;
    std::list<Point> randPoints;
    ball B0;
    Point q(n);
    PolyBall zb_it;
    //get_first_ball(Polytope &P, ball &B0, NT &ratio, NT &radius, const NT &lb, const NT &ub, const NT &alpha,NT &rmax)
    if( !get_first_ball<RNGType>(P, B0, ratio, radius, lb, ub, alpha, rmax) ) {
        return false;
    }
    //NT rad_min = B0.radius();
    //std::cout<<"first ball computed"<<std::endl;
    ratio0 = ratio;
    rand_point_generator(P, q, Ntot, var.walk_steps, randPoints, var);
    //std::cout<<Ntot<<" points sampled from P"<<std::endl;

    if (check_convergence<Point>(B0, randPoints, lb, ub, fail, ratio, nu, alpha, false, true)) {
        ratios.push_back(ratio);
        BallSet.push_back(B0);
        ratios.push_back(ratio0);
        //if(print) std::cout<<"one ball and ratio = "<<ratio<<std::endl;
        return true;
    }
    //if(print) std::cout<<"not the last ball, ratio = "<<ratio<<std::endl;
    if ( !get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, alpha, nu) ){
        return false;
    }
    //if(print) std::cout<<"number of balls = "<<BallSet.size()+1<<std::endl;

    while (true) {
        zb_it = PolyBall(P, BallSet[BallSet.size()-1]);
        q=Point(n);
        randPoints.clear();
        //zb_it.comp_diam(var.diameter);
        rand_point_generator(zb_it, q, Ntot, var.walk_steps, randPoints,var);
        //std::cout<<"N points sampled from BP"<<std::endl;

        if (check_convergence<Point>(B0, randPoints, lb, ub, fail, ratio, nu, alpha, false, true)) {
            ratios.push_back(ratio);
            BallSet.push_back(B0);
            ratios.push_back(ratio0);
            //if(print) std::cout<<"annealing stoped, last ratio = "<<ratio<<std::endl;
            return true;
        }
        //if(print) std::cout<<"compute new ball"<<std::endl;
        if ( !get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, alpha, nu) ) {
            return false;
        }
        //if(print) std::cout<<"number of balls = "<<BallSet.size()+1<<std::endl;
    }
}


#endif

