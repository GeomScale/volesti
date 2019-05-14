// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_ANNEALING_H
#define BALL_ANNEALING_H




template <class Hpolytope, class PointList, typename NT>
void check_converg22(Hpolytope &P, PointList &randPoints, NT p_test, bool &done, bool &too_few, NT &ratio, NT up_lim, int nu, bool print) {

    typedef typename Hpolytope::PolytopePoint Point;
    std::vector<NT> ratios;
    NT countsIn = 0.0;
    int m = randPoints.size()/nu;

    int i = 1;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){
        if (P.is_in(*pit)==-1) {
            countsIn += 1.0;
        }
        if (i % m == 0) {
            //if (print) std::cout<<"ratio = "<<countsIn/120.0<<std::endl;
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
        }
        ratio = p_mval; // for test only
    } else {
        too_few = true;
        ratio = p_mval; // for test only
    }

}


template <class PointList, class HPolytope, typename NT, class Parameter>
bool is_last_zonoball22(PointList randPoints, HPolytope &HP, NT &ratio, NT p_value, int nu, Parameter var){

    NT countIn = 0.0;
    std::vector<NT> ratios;
    int m = randPoints.size()/nu;
    bool print = var.verbose;

    typename PointList::iterator rpit = randPoints.begin();
    int i=1;
    for ( ;  rpit!=randPoints.end(); ++rpit, ++i) {

        if (HP.is_in(*rpit)==-1) {
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
    //if (print) std::cout<<"check for last hpoly, p_mval = "<<p_mval<<std::endl;
    if (p_mval > p_value + t_value*(p_varval/std::sqrt(NT(ni)))) {
        ratio = p_mval;
        return true;
    }
    return false;

}

template <class Zonotope, class HPolytope, class VT, class PointList, typename NT, class Parameters>
void get_next_zonoball22(Zonotope &Z, std::vector<HPolytope> &HPolySet,
                         HPolytope &HP2, VT Zs_max, VT Zs_min, PointList randPoints,
                        std::vector<NT> &ratios, NT p_value, NT up_lim, int nu, Parameters &var){

    //typename typedef ZonoBall::ball ball;
    typedef typename Zonotope::PolytopePoint Point;
    int n = var.n;
    bool done, too_few;
    bool print = var.verbose;

    NT rad2=0.0;
    NT rad1=0.0, rad;
    NT pnorm, ratio;

    VT Zmed(Zs_max.size());
    NT med, u = 1.0, l = 0.0;

    while (true) {
        med = (u + l) * 0.5;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HP2.set_vec(Zmed);
        done = false;
        too_few = false;

        check_converg22(HP2, randPoints, p_value, done, too_few, ratio, up_lim, nu, false);

        if(done){
            HPolySet.push_back(HP2);
            ratios.push_back(ratio);
            return;
        }

        if(print) std::cout<<"med = "<<med<<" ratio = "<<ratio<<std::endl;
        if(too_few) {
            l = med;
        } else {
            u = med;
        }

    }

}

template <class ZonoHP,class Zonotope, class HPolytope, class VT, class Parameters, typename NT>
void get_sequence_of_zonopolys(Zonotope &Z, HPolytope &HP, std::vector<HPolytope> &HPolySet,
                               VT Zs_max, std::vector<NT> &ratios, int Ntot, int nu,
                               NT &p_value, NT up_lim, Parameters &var) {

    bool print = var.verbose;
    typedef typename Zonotope::PolytopePoint Point;
    typedef typename Zonotope::MT MT;
    //MT Q0 = Z.get_Q0().transpose();
    //MT G = Z.get_mat().transpose();
    int n = var.n;
    //int k = Z.num_of_generators();
    //bool done = false;
    NT ratio;
    std::list<Point> randPoints;
    Point q(n);
    //HnRsteps = 0.0;
    //std::cout<<"Ntot = "<<Ntot<<" nu = "<<nu<<std::endl;
    rand_point_generator(Z, q, Ntot, var.walk_steps, randPoints, var);
    //HnRsteps += NT(Ntot)*NT(var.walk_steps);
    //PointSets.push_back(randPoints);
    HPolytope HP2 = HP;
    if (is_last_zonoball22(randPoints, HP, ratio, p_value, nu, var)) {
        ratios.push_back(ratio);
        if(print) std::cout<<"last hpoly and ratio = "<<ratio<<std::endl;
        return;
    }
    if(print) std::cout<<"not the last hpoly"<<std::endl;
    get_next_zonoball22(Z, HPolySet, HP2, Zs_max, HP.get_vec(), randPoints, ratios, p_value, up_lim, nu, var);
    if(print) std::cout<<"get first hpoly"<<std::endl;

    ZonoHP ZHP2;
    VT Zs_min = HP.get_vec();

    while (true) {
        //HP2 = HPolySet[HPolySet.size()-1];
        ZHP2 = ZonoHP(Z,HP2);
        q=Point(n);
        randPoints.clear();
        rand_point_generator(ZHP2, q, Ntot, var.walk_steps, randPoints, var);
        //HnRsteps += NT(Ntot)*NT(var.walk_steps);
        //PointSets.push_back(randPoints);
        if (is_last_zonoball22(randPoints, HP, ratio, p_value, nu, var)) {
            ratios.push_back(ratio);
            if(print) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
            return;
        }
        get_next_zonoball22(Z, HPolySet, HP2, HP2.get_vec(), Zs_min, randPoints, ratios, p_value, up_lim, nu, var);
        if(print) std::cout<<"get hpoly"<<std::endl;
    }


}

#endif
