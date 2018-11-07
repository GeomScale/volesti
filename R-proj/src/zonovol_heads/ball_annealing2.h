// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALINGTWO_H
#define BALL_ANNEALINGTWO_H


template <typename NT>
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
}

template <class Point, class ball, class PointList, typename NT>
void check_converg2(ball &P, PointList &randPoints, NT p_test, bool &done, bool &too_few, NT &ratio, NT up_lim, bool print) {

    //typedef typename ball::BallPoint Point;
    std::vector<NT> ratios;
    NT countsIn = 0.0;

    int i = 1;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (P.is_in(*pit)==-1) {
            countsIn += 1.0;
        }
        if (i % 120 == 0) {
            //if (print) std::cout<<"ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countsIn/120.0);
            countsIn = 0.0;
        }
    }

    std::pair<NT,NT> mv = getMeanVariance(ratios);
    NT t_value = 0.700;
    NT p_mval = mv.first;
    NT p_varval = mv.second;
    int ni = ratios.size();
    //NT p_test = a;

    //if (print) std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
    if (p_mval > p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))) {
        if (p_mval < (up_lim) + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))) {
            done= true;
            ratio = p_mval;
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
bool is_last_zonoball(PointList randPoints, ball &B0, NT  &ratio, Parameter var){

    NT countIn = 0.0;
    std::vector<NT> ratios;


    typename PointList::iterator rpit = randPoints.begin();
    int i=1;
    for ( ;  rpit!=randPoints.end(); ++rpit, ++i) {
        if(B0.is_in(*rpit)==-1) {
            countIn = countIn + 1.0;
        }
        if (i % 120 == 0) {
            //if (print) std::cout<<"ratio = "<<countsIn/120.0<<std::endl;
            ratios.push_back(countIn/120.0);
            countIn = 0.0;
        }
    }

    std::pair<NT,NT> mv = getMeanVariance(ratios);
    NT t_value = 0.700;
    NT p_mval = mv.first;
    NT p_varval = mv.second;
    int ni = ratios.size();
    //NT p_test = a;

    //if (print) std::cout<<"mean must be greater than = "<<p_test + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
    if (p_mval > 0.1 + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))) {
        ratio = p_mval;
        return true;
    }
    return false;

}

template <class ball, class Zonotope, class PointList, typename NT, class Parameters>
void get_next_zonoball(Zonotope &Z, std::vector<ball> &BallSet,
                       PointList &randPoints, std::vector<NT> &ratios, Parameters &var){

    //typename typedef ZonoBall::ball ball;
    typedef typename Zonotope::PolytopePoint Point;
    int n = var.n;
    bool done, too_few;

    NT rad2=0.0;
    NT rad1=0.0, rad;
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

        check_converg2<Point>(Biter, randPoints, 0.1, done, too_few, ratio, 0.15, false);

        if(done){
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return;
        }

        if(too_few) {
            rad1 = rad;
        } else {
            rad2 = rad;
        }
    }

}

template <class RNGType,class ball, class Zonotope, class Parameters, typename NT>
void get_first_ball(Zonotope &Z, ball &B0, NT &ratio, NT radius, Parameters &var){

    typedef typename Zonotope::PolytopePoint Point;
    NT rad2 = 10*radius;
    NT rad1 = radius;
    int n = var.n;
    bool done, too_few;

    ball BallIter;
    Point p(n);
    std::list<Point> randPoints;
    Point center(n);

    while(true) {

        randPoints.clear();
        for (int i = 0; i < 1200; ++i) {
            p = get_point_in_Dsphere<RNGType, Point>(n, rad2);
            randPoints.push_back(p);
        }
        done = false; too_few = false;
        check_converg2<Point>(Z, randPoints, 0.1, done, too_few, ratio, 0.15, false);
        std::cout<<"rad2 = "<<rad2<<std::endl;
        std::cout<<"ratio = "<<ratio<<std::endl;

        if(done) {
            BallIter = ball(center, rad2*rad2);
            //zb = ZonoBall(Z,BallIter);
            B0 = BallIter;
            return;
        }

        if (too_few) {
            break;
        }
        rad1 = rad2;
        rad2 = rad2 + 10*radius;
    }

    NT rad_med;

    while(true) {
        rad_med = 0.5*(rad1+rad2);
        randPoints.clear();
        for (int i = 0; i < 1200; ++i) {
            p = get_point_in_Dsphere<RNGType, Point>(n, rad_med);
            randPoints.push_back(p);
        }
        done = false; too_few = false;
        check_converg2<Point>(Z, randPoints, 0.1, done, too_few, ratio, 0.15, false);

        std::cout<<"rad_med = "<<rad_med<<std::endl;
        std::cout<<"ratio = "<<ratio<<std::endl;
        if(done) {
            BallIter = ball(center, rad_med*rad_med);
            //zb = ZonoBall(Z,BallIter);
            B0 = BallIter;
            return;
        }

        if (too_few) {
            rad2 = rad_med;
        } else {
            rad1 = rad_med;
        }

    }

}

template <class ZonoBall, class RNGType,class ball, class Zonotope, class PointList, class Parameters, typename NT>
void get_sequence_of_zonoballs(Zonotope &Z, std::vector<ball> &BallSet, ball &B0, NT &ratio0,
                               std::vector<PointList> &PointSets, std::vector<NT> &ratios,
                               NT &p_value, NT radius, Parameters &var) {


    typedef typename Zonotope::PolytopePoint Point;
    typedef typename Zonotope::MT MT;
    MT Q0 = Z.get_Q0().transpose();
    MT G = Z.get_mat().transpose();
    int n = var.n;
    //int k = Z.num_of_generators();
    //bool done = false;
    NT ratio;
    PointList randPoints;
    Point q(n);
    ZonoBall zb_it;
    //ball B0;
    get_first_ball<RNGType>(Z, B0, ratio, radius, var);
    ratio0 = ratio;

    rand_point_generator(Z, q, 1200, 1, randPoints, var);
    PointSets.push_back(randPoints);
    if (is_last_zonoball(randPoints, B0, ratio, var)) {
        ratios.push_back(ratio);
        std::cout<<"one ball and ratio = "<<ratio<<std::endl;
        return;
    }
    std::cout<<"not the last ball, ratio = "<<ratio<<std::endl;
    get_next_zonoball(Z, BallSet, randPoints, ratios, var);

    ball Biter;

    while (true) {
        Biter = BallSet[BallSet.size()-1];
        zb_it = ZonoBall(Z,Biter);
        q=Point(n);
        randPoints.clear();
        rand_point_generator(zb_it, q, 1200, 1, randPoints, var);
        PointSets.push_back(randPoints);
        if (is_last_zonoball(randPoints, B0, ratio, var)) {
            ratios.push_back(ratio);
            //std::cout<<"number of balls = "<<ZonoBallSet.size()<<std::endl;
            return;
        }
        get_next_zonoball(Z, BallSet, randPoints, ratios, var);
    }


}

#endif
