// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef BALL_HPOLY_ANNEALING_H
#define BALL_HPOLY_ANNEALING_H


template <class VPolytope, class BallPoly, typename NT, class Parameters>
void get_first_ballpoly(VPolytope &VP, BallPoly &BP, NT lb, NT &up_lim, NT &ratio, Parameters &var){

    typedef typename VPolytope::PolytopePoint Point;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;
    //MT G = P.get_mat().transpose();
    //MT A = HP.get_mat();
    //int kk = G.cols();
    VT Zs_max = BP.get_vec();
    VT Zs_min = VT::Zero(BP.num_of_hyperplanes());
    BallPoly BPiter=BP;

    int n = VP.dimension(), m = Zs_min.size();
    int N = 1200;
    Point q(n);
    bool too_few, print = false;
    std::list<Point> randPoints;

    NT l=0.0, u=1.0, med;
    VT  Zmed = Zs_max;
    int count =0;
    Parameters variter = var;
    NT a, rad_med = BP.radius(), r, a_initial = Zmed(0);
    while(true) {

        count++;
        q=Point(n);
        med = (u + l) * 0.5;
        std::cout<<"u = "<<u<<" l = "<<l<<" med = "<<med<<std::endl;
        a = Zmed(0);
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        BPiter.set_vec(Zmed);
        r = Zmed(0)/a;
        //r = 1.0/r;
        BPiter.set_radius(BPiter.radius()*r);
        //variter.che_rad = HPiter.ComputeInnerBall().second;

        randPoints.clear();
        rand_point_generator(BPiter, q, 1200, 10+2*n, randPoints, variter);
        var.MemLps = var.MemLps + 1200.0;
        too_few = false;

        if(check_converg001<Point>(VP, randPoints, lb, up_lim, too_few, ratio, 10, 0.2, true, false)) {
            BP.set_vec(Zmed);
            r = Zmed(0)/a_initial;
            //r = 1.0/r;
            BP.set_radius(BPiter.radius()*r);
            return;
        }

        if (too_few) {
            u = med;
        } else {
            l = med;
        }
        if(med>0.9) {
            BP.set_vec(Zmed);
            r = Zmed(0)/a_initial;
            //r = 1.0/r;
            BP.set_radius(BPiter.radius()*r);
            return;
        }
        if(u-l<0.0000000001) {
            std::cout << "fail to find first hpoly... repeat proccess" << std::endl;
            //std::cout<<"origin is in = "<<P.is_in(Point(n))<<std::endl;
            u=1.0;
            l=0.0;
        }
    }
}


template <class Zonotope, class BallPoly, class VT, class PointList, typename NT>
void get_next_ballpoly(Zonotope &Z, std::vector<BallPoly> &BPolySet,
                       BallPoly &BP2, VT Zs_max, VT Zs_min, PointList randPoints,
                         std::vector<NT> &ratios, NT p_value, NT up_lim, int nu, NT alpha, bool first){

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
    bool too_few;
    bool print = false;

    NT rad2=0.0;
    NT rad1=0.0, rad;
    NT pnorm, ratio;

    VT Zmed(Zs_max.size());
    if (first) {
        Zmed = Zs_min;
    } else {
        Zmed = Zs_max;
    }
    NT med, u = 1.0, l = 0.0, a,r, a_initial = Zmed(0);

    while (true) {
        med = (u + l) * 0.5;
        a = Zmed(0);
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        BP2.set_vec(Zmed);
        r = Zmed(0)/a;
        //r = 1.0/r;
        BP2.set_radius(BP2.radius()*r);
        too_few = false;

        if(check_converg001<Point>(BP2, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, false)){
            BPolySet.push_back(BP2);
            ratios.push_back(ratio);
            return;
        }
        if(too_few) {
            l = med;
        } else {
            u = med;
        }
    }
}


template <class ZonoHP,class Vpolytope, class BallPoly, class Parameters, typename NT>
void get_sequence_of_vpoly_ballpolys(Vpolytope &VP, BallPoly &BP, std::vector<BallPoly> &BPolySet,
                                  std::vector<NT> &ratios, int Ntot, int nu,
                                  NT &p_value, NT up_lim, NT alpha, Parameters &var, Parameters &var2) {

    bool print = var.verbose, too_few=false;
    typedef typename Vpolytope::PolytopePoint Point;
    typedef typename Vpolytope::MT MT;
    typedef typename Vpolytope::VT VT;

    int n = var.n;
    VT Zs_max = VT::Ones(BP.num_of_hyperplanes());

    NT ratio;
    std::list<Point> randPoints;
    Point q(n);
    std::cout<<"sample = "<<Ntot<<" points"<<std::endl;
    //std::cout<<"walk_step = "<<var.walk_steps<<std::endl;
    rand_point_generator(VP, q, Ntot, var.walk_steps, randPoints, var);
    var.TotSteps = var.TotSteps + NT(Ntot);
    std::cout<<"points sampled"<<std::endl;
    if (check_converg001<Point>(BP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true)) {
        ratios.push_back(ratio);
        if(print) std::cout<<"last hpoly and ratio = "<<ratio<<std::endl;
        return;
    }
    BallPoly BP2 = BP;
    if(print) std::cout<<"not the last hpoly"<<std::endl;
    get_next_ballpoly(VP, BPolySet, BP2, Zs_max, BP.get_vec(), randPoints, ratios, p_value, up_lim, nu, alpha, true);
    if(print) std::cout<<"get first hpoly"<<std::endl;

    ZonoHP ZHP2;
    VT Zs_min = BP.get_vec();

    while (true) {

        ZHP2 = ZonoHP(VP,BP2);
        q=Point(n);
        randPoints.clear();

        var.diameter = 2.0*BP2.radius();
        std::cout<<"[annealing] diameter = "<<var.diameter<<std::endl;
        //var.diameter = diams_inter[diams_inter.size()-1];
        std::cout<<"sampling N points from ZHP2, walk_length ="<<var.walk_steps<<std::endl;
        rand_point_generator(ZHP2, q, Ntot, var.walk_steps, randPoints, var);
        var.TotSteps = var.TotSteps + NT(Ntot);
        std::cout<<"Checking convergence"<<std::endl;
        if (check_converg001<Point>(BP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true)) {
            ratios.push_back(ratio);
            if(print) std::cout<<"number of hpolys = "<<BPolySet.size()<<std::endl;
            return;
        }
        std::cout<<"getting a new hpoly"<<std::endl;
        get_next_ballpoly(VP, BPolySet, BP2, BP2.get_vec(), Zs_min, randPoints, ratios, p_value, up_lim, nu, alpha, false);
        if(print) std::cout<<"got hpoly"<<std::endl;
    }


}


#endif
