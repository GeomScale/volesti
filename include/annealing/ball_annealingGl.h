// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALINGGL_H
#define BALL_ANNEALINGGL_H

/*template <typename NT>
std::pair<NT, NT> getMeanVariance(std::vector<NT>& vec) {
    NT mean = 0, M2 = 0, variance = 0, delta;
    typedef typename std::vector<NT>::iterator viterator;

    unsigned int i=0;
    viterator vecit = vec.begin();
    for( ; vecit!=vec.end(); vecit++, i++){
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }

    return std::pair<NT, NT> (mean, variance);
}*/

template <class Point, class ball, class PointList, typename NT>
void check_converg001(ball &P, PointList &randPoints, NT p_test, bool &done, bool &too_few, NT &ratio,
                      NT up_lim, int nu, bool print, std::vector<NT> &Ci_ratios) {

    //typedef typename ball::BallPoint Point;
    std::vector<NT> ratios;
    NT countsIn = 0.0;
    int m = randPoints.size()/nu;
    NT pr = 0.99, rm , rs;
    NT p = (1.0+pr)/2.0;
    NT zp = std::sqrt(2.0)*boost::math::erf_inv(2.0*p - 1.0);
    std::pair<NT,NT> mv;

    int i = 1, count_sets=0;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (P.is_in(*pit)==-1) {
            countsIn += 1.0;
        }
        if (i % m == 0) {
            //count_sets++;
            if (print) std::cout<<"ratio = "<<countsIn/m<<std::endl;
            ratios.push_back(countsIn/m);
            countsIn = 0.0;
            mv = getMeanVariance(ratios);
            rm = mv.first; rs = mv.second;
            //std::cout<<"rm = "<<rm<<"rs = "<<rs<<"rm+zp*rs = "<<rm+zp*rs<<"rm-zp*rs = "<<rm-zp*rs<<std::endl;
            if (rm+zp*rs<(p_test-0.01)) {
                too_few = true;
                return;
            } else if (rm-zp*rs>(up_lim+0.02)) {
                return;
            }
        }

        //if (count_sets==1) {
        //if(ratios[0]<0.04) {
        //too_few=true;
        //return;
        //} else if(ratios[0]>0.5) {
        //return;
        //}
        //}
    }

    mv = getMeanVariance(ratios);
    NT t_value = 0.700;
    NT p_mval = mv.first;
    NT p_varval = std::sqrt(mv.second);
    int ni = ratios.size();
    //NT p_test = a;

    //if (print) std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
    if (p_mval > p_test + t_value*(p_varval/std::sqrt(NT(ni)))) {
        if (p_mval < (up_lim) + t_value*(p_varval/std::sqrt(NT(ni)))) {
            done= true;
            ratio = p_mval;
            Ci_ratios = ratios;
            //std::cout<<"ni = "<<ni<<std::endl;
            //std::cout<<"ratio done = "<<ratio<<" var = "<<p_varval<<" "<<(up_lim) + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
        }
        ratio = p_mval; // for test only
    } else {
        too_few = true;
        ratio = p_mval; // for test only
    }

}

template <class Point, class ball, class PointList, typename NT>
void check_converg2(ball &P, PointList &randPoints, NT p_test, bool &done, bool &too_few, NT &ratio, NT up_lim,
                    int nu, bool print, std::vector<NT> &Ci_ratios) {

    //typedef typename ball::BallPoint Point;
    std::vector<NT> ratios;
    NT countsIn = 0.0;
    int m = randPoints.size()/nu;

    int i = 1;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (P.is_in(*pit)==-1) {
            countsIn += 1.0;
        }
        if (i % m == 0) {
            if (print) std::cout<<"ratio = "<<countsIn/m<<std::endl;
            ratios.push_back(countsIn/m);
            countsIn = 0.0;
        }
    }

    std::pair<NT,NT> mv = getMeanVariance(ratios);
    NT t_value = 0.700;
    NT p_mval = mv.first;
    NT p_varval = std::sqrt(mv.second);
    int ni = ratios.size();
    //NT p_test = a;

    //if (print) std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
    if (p_mval > p_test + t_value*(p_varval/std::sqrt(NT(ni)))) {
        if (p_mval < (up_lim) + t_value*(p_varval/std::sqrt(NT(ni)))) {
            done= true;
            ratio = p_mval;
            Ci_ratios = ratios;
            //std::cout<<"ni = "<<ni<<std::endl;
            //std::cout<<"ratio done = "<<ratio<<" var = "<<p_varval<<" "<<(up_lim) + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
        }
        ratio = p_mval; // for test only
    } else {
        too_few = true;
        ratio = p_mval; // for test only
    }

}


template <class PointList, class ball, typename NT, class Parameter>
bool is_last_zonoball(PointList randPoints, ball &B0, NT  &ratio, NT p_test, int nu, Parameter var, std::vector<NT> &Ci_ratios){

    NT countIn = 0.0;
    std::vector<NT> ratios;
    int m = randPoints.size()/nu;

    //std::cout<<"helloo"<<std::endl;
    typename PointList::iterator rpit = randPoints.begin();
    int i=1;
    for ( ;  rpit!=randPoints.end(); ++rpit, ++i) {
        if(B0.is_in(*rpit)==-1) {
            countIn = countIn + 1.0;
        }
        if (i % m == 0) {
            //if (print) std::cout<<"ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countIn/m);
            countIn = 0.0;
        }
    }

    std::pair<NT,NT> mv = getMeanVariance(ratios);
    NT t_value = 0.700;
    NT p_mval = mv.first;
    NT p_varval = std::sqrt(mv.second);
    int ni = ratios.size();
    //NT p_test = a;

    //if (print) std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
    if (p_mval > p_test + t_value*(p_varval/std::sqrt(NT(ni)))) {
        ratio = p_mval;
        Ci_ratios = ratios;
        //std::cout<<"last zonoball"<<std::endl;
        return true;
    }
    //std::cout<<"not the last zonoball"<<std::endl;
    return false;

}

template <class ball, class Zonotope, class PointList, typename NT, class Parameters>
std::vector<NT> get_next_zonoball(Zonotope &Z, std::vector<ball> &BallSet,
                       PointList &randPoints, NT rad_min, std::vector<NT> &ratios, NT &p_value, NT up, int nu,
                       Parameters &var){

    //typename typedef ZonoBall::ball ball;
    typedef typename Zonotope::PolytopePoint Point;
    std::vector<NT> Ci_ratios;
    int n = var.n;
    bool done, too_few;
    bool print = var.verbose;

    NT rad2=0.0;
    NT rad1=rad_min, rad;
    NT pnorm, ratio;
    Point center(n);

    typename PointList::iterator rpit = randPoints.begin();
    for ( ;  rpit!=randPoints.end(); ++rpit) {
        pnorm = (*rpit).squared_length();
        if(pnorm >rad2){
            rad2=pnorm;
        }
    }
    ball Biter;
    rad2=std::sqrt(rad2);

    while (true) {
        rad = 0.5*(rad1+rad2);
        Biter = ball(center, rad*rad);
        done = false;
        too_few = false;

        check_converg2<Point>(Biter, randPoints, p_value, done, too_few, ratio, up, nu, false, Ci_ratios);

        //if(print) std::cout<<"rad = "<<rad<<" ratio = "<<ratio<<std::endl;
        if(done){
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return Ci_ratios;
        }

        if(too_few) {
            rad1 = rad;
        } else {
            rad2 = rad;
        }
    }

}

template <class RNGType,class ball, class Zonotope, class Parameters, typename NT>
std::vector<NT> get_first_ball(Zonotope &Z, ball &B0, NT &ratio, NT radius, NT lb, NT up, Parameters &var, NT rmax){

    typedef typename Zonotope::PolytopePoint Point;
    std::vector<NT> Ci_ratios;
    int n = var.n;
    NT rad2;
    bool bis_int = false;
    bool print = var.verbose;
    if(rmax>0.0) {
        rad2 = rmax;
        std::list<Point> randPoints2;
        randPoints2.clear();
        Point pp(n);
        for (int i = 0; i < 1200; ++i) {
            pp = get_point_in_Dsphere<RNGType, Point>(n, rad2);
            randPoints2.push_back(pp);
        }
        //ballsteps +=1200.0;
        bool done2 = false, too_few2 = false;
        check_converg001<Point>(Z, randPoints2, lb, done2, too_few2, ratio, up, 10, false, Ci_ratios);
        if (done2 || !too_few2) {
            B0 = ball(Point(n), rad2*rad2);
            if(print) std::cout<<"rmax is enclosing and ok for rejection.., rmax = "<<rmax<<" rmin = "<<radius<<std::endl;
            return Ci_ratios;
        } else {
            if(print) std::cout<<"rmax is NOT OK for rejection..rmax = "<<rmax<<" rmin = "<<radius<<std::endl;
        }
        bis_int = true;
    } else {
        rad2 = 2 * std::sqrt(NT(n)) * radius;
    }
    NT rad1 = radius;

    bool done, too_few;

    ball BallIter;
    Point p(n);
    std::list<Point> randPoints;
    Point center(n);

    while(!bis_int) {

        randPoints.clear();
        for (int i = 0; i < 1200; ++i) {
            p = get_point_in_Dsphere<RNGType, Point>(n, rad2);
            randPoints.push_back(p);
        }
        //ballsteps +=1200.0;
        done = false; too_few = false;
        check_converg001<Point>(Z, randPoints, lb, done, too_few, ratio, up, 10, false, Ci_ratios);
        if(print) std::cout<<"rad2 = "<<rad2<<std::endl;
        if(print) std::cout<<"ratio = "<<ratio<<std::endl;

        if(done) {
            BallIter = ball(center, rad2*rad2);
            //zb = ZonoBall(Z,BallIter);
            B0 = BallIter;
            return Ci_ratios;
        }

        if (too_few) {
            break;
        }
        rad1 = rad2;
        rad2 = rad2 + 2*std::sqrt(NT(n))*radius;
    }

    NT rad_med;

    while(true) {
        rad_med = 0.5*(rad1+rad2);
        randPoints.clear();
        for (int i = 0; i < 1200; ++i) {
            p = get_point_in_Dsphere<RNGType, Point>(n, rad_med);
            randPoints.push_back(p);
        }
        //ballsteps +=1200.0;
        done = false; too_few = false;
        check_converg001<Point>(Z, randPoints, lb, done, too_few, ratio, up, 10, false, Ci_ratios);

        if(print) std::cout<<"rad_med = "<<rad_med<<std::endl;
        if(print) std::cout<<"ratio = "<<ratio<<std::endl;
        if(done) {
            if(print) std::cout<<"rad_med = "<<rad_med<<std::endl;
            BallIter = ball(center, rad_med*rad_med);
            //zb = ZonoBall(Z,BallIter);
            B0 = BallIter;
            return Ci_ratios;
        }

        if (too_few) {
            rad2 = rad_med;
        } else {
            rad1 = rad_med;
        }

    }

}

template <class ZonoBall, class RNGType,class ball, class Zonotope, class Parameters, typename NT>
void get_sequence_of_zonoballs(Zonotope &Z, std::vector<ball> &BallSet, ball &B0, NT &ratio0,
                               std::vector<NT> &ratios, int Ntot, int nu,
                               NT &p_value, NT up, NT radius, Parameters &var,
                               NT B0_radius = 0.0, NT rmax = 0.0) {


    typedef typename Zonotope::PolytopePoint Point;
    typedef typename Zonotope::MT MT;
    bool print = var.verbose;
    //MT Q0 = Z.get_Q0().transpose();
    //MT G = Z.get_mat().transpose();
    //HnRSteps = 0.0;
    int n = var.n;
    //int k = Z.num_of_generators();
    //bool done = false;
    NT ratio;
    std::list<Point> randPoints;
    Point q(n);
    ZonoBall zb_it;
    //ball B0;
    //ΝΤ ballSteps = 0.0;
    std::vector<NT> Ci_ratios;
    if (B0_radius==0.0) {
        Ci_ratios = get_first_ball<RNGType>(Z, B0, ratio, radius, p_value, up, var, rmax);
        //all_ratios.push_back(Ci_ratios);
        ratio0 = ratio;
    } else {
        B0 = ball(Point(n), B0_radius*B0_radius);
    }

    //std::cout<<"Ntot = "<<Ntot<<" nu = "<<nu<<std::endl;
    rand_point_generator(Z, q, Ntot, var.walk_steps, randPoints, var);
    //std::cout<<"points sampled"<<std::endl;
    //HnRSteps += NT(Ntot)*NT(var.walk_steps);
    //PointSets.push_back(randPoints);
    if (is_last_zonoball(randPoints, B0, ratio, p_value, nu, var, Ci_ratios)) {
        //all_ratios.push_back(Ci_ratios);
        ratios.push_back(ratio);
        if(print) std::cout<<"one ball and ratio = "<<ratio<<std::endl;
        return;
    }
    if(print) std::cout<<"not the last ball, ratio = "<<ratio<<std::endl;
    Ci_ratios = get_next_zonoball(Z, BallSet, randPoints, B0.radius(), ratios, p_value, up, nu, var);
    //all_ratios.push_back(Ci_ratios);
    if(print) std::cout<<"number of balls = "<<BallSet.size()<<std::endl;

    ball Biter;

    while (true) {
        Biter = BallSet[BallSet.size()-1];
        zb_it = ZonoBall(Z,Biter);
        q=Point(n);
        randPoints.clear();
        rand_point_generator(zb_it, q, Ntot, var.walk_steps, randPoints, var);
        //HnRSteps += NT(Ntot)*NT(var.walk_steps);
        //PointSets.push_back(randPoints);
        if (is_last_zonoball(randPoints, B0, ratio, p_value, nu, var, Ci_ratios)) {
            //all_ratios.push_back(Ci_ratios);
            ratios.push_back(ratio);
            //std::cout<<"number of balls = "<<ZonoBallSet.size()<<std::endl;
            return;
        }
        Ci_ratios = get_next_zonoball(Z, BallSet, randPoints, B0.radius(), ratios, p_value, up, nu, var);
        //all_ratios.push_back(Ci_ratios);
        if(print) std::cout<<"number of balls = "<<BallSet.size()<<std::endl;
    }


}


#endif
