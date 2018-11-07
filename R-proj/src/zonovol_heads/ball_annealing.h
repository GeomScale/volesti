// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALING_H
#define BALL_ANNEALING_H

template <class ball, class PointList, typename NT>
void check_converg2(ball &P, PointList &randPoints, NT p_test, bool &done, bool &too_few, NT &ratio, NT up_lim, bool print) {

    typedef typename ball::BallPoint Point;
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
        }
        ratio = p_mval; // for test only
    } else {
        too_few = true;
        ratio = p_mval; // for test only
    }

}


template <class PointList, class HPolytope, typename NT, class Parameter>
bool is_last_zonoball(PointList randPoints, HPolytope &HP, NT &ratio, Parameter var){

    NT countIn = 0.0;
    std::vector<NT> ratios;


    typename PointList::iterator rpit = randPoints.begin();
    int i=1;
    for ( ;  rpit!=randPoints.end(); ++rpit, ++i) {

        if (HP.is_in(*rpit)==-1) {
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

template <class ball, class Zonotope, class ZonoBall, class PointList, typename NT, class Parameters>
void get_next_zonoball(Zonotope &Z, std::vector<ZonoBall> &ZonoBallSet,
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

        check_converg2(Biter, randPoints, 0.1, done, too_few, ratio, 0.15, false);

        if(done){
            ZonoBallSet.push_back(ZonoBall(Z, Biter));
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

template <class ball, class Zonotope, class HPolytope, class ZonoBall, class PointList, class Parameters, typename NT>
void get_sequence_of_zonoballs(Zonotope &Z, HPolytope &HP, std::vector<ZonoBall> &ZonoBallSet,
                               std::vector<PointList> &PointSets, std::vector<NT> &ratios,
                               NT &p_value, Parameters &var, NT &delta, std::vector<NT> &Zs, bool relaxed) {


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
    rand_point_generator(Z, q, 1200, 1, randPoints, var);
    PointSets.push_back(randPoints);
    if (is_last_zonoball(randPoints, HP, ratio, var)) {
        ratios.push_back(ratio);
        std::cout<<"last ball and ratio = "<<ratio<<std::endl;
        return;
    }
    std::cout<<"not the last ball"<<std::endl;
    get_next_zonoball<ball>(Z, ZonoBallSet, randPoints, ratios, var);

    ZonoBall ZBiter;

    while (true) {
        ZBiter = ZonoBallSet[ZonoBallSet.size()-1];
        q=Point(n);
        randPoints.clear();
        rand_point_generator(ZBiter, q, 1200, 1, randPoints, var);
        PointSets.push_back(randPoints);
        if (is_last_zonoball(randPoints, HP, ratio, var)) {
            ratios.push_back(ratio);
            std::cout<<"number of balls = "<<ZonoBallSet.size()<<std::endl;
            return;
        }
        get_next_zonoball<ball>(Z, ZonoBallSet, randPoints, ratios, var);
    }


}

#endif
