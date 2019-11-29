// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef HPOLY_ANN_PROJ_H
#define HPOLY_ANN_PROJ_H




template <class VPolytope, class HPolytope, typename NT, class Parameters>
void get_first_poly33(VPolytope &VP, HPolytope &HP, NT lb, NT &up_lim, NT &ratio, Parameters &var){

    typedef typename VPolytope::PolytopePoint Point;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;
    //MT G = P.get_mat().transpose();
    //MT A = HP.get_mat();
    //int kk = G.cols();
    VT Zs_max = HP.get_vec();
    VT Zs_min = VT::Zero(HP.num_of_hyperplanes());
    HPolytope HPiter=HP;

    int n = VP.dimension(), m = Zs_min.size();
    int N = 1200,count;
    Point q(n);
    bool too_few, print = false;
    std::list<Point> randPoints;

    NT l=0.0, u=1.0, med;
    VT  Zmed(m);
    int count2 =0;
    Parameters variter = var;
    while(true) {

        count++;
        q=Point(n);
        med = (u + l) * 0.5;
        std::cout<<"u = "<<u<<" l = "<<l<<" med = "<<med<<std::endl;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HPiter.set_vec(Zmed);
        //variter.che_rad = HPiter.ComputeInnerBall().second;

        randPoints.clear();
        rand_point_generator(HPiter, q, 1200, n, randPoints, variter);
        var.MemLps = var.MemLps + 1200.0;
        too_few = false;

        check_converg001<Point>(VP, randPoints, lb, up_lim, too_few, ratio, 10, 0.2, false, false);
        //count2 = 0;
        //for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit) {
            //if(VP.is_in(*pit)==-1) count2++;
        //}
        //ratio = NT(count2)/NT(1200.0);
        std::cout << "ratio = "<<ratio<< std::endl;
        if(ratio>0.09 && ratio<0.12) {
            HP.set_vec(Zmed);
            return;
        } else if(ratio>0.12){
            l=med;
        } else {
            u=med;
        }
        std::cout << "ratio = "<<ratio<< std::endl;

        //if (too_few) {
            //u = med;
        //} else {
            //l = med;
        //}
        if(med>0.9) {
            HP.set_vec(Zmed);
            return;
        }
        if(u-l<0.00000000000001) {
            std::cout << "fail to find first hpoly... repeat proccess" << std::endl;
            //std::cout<<"origin is in = "<<P.is_in(Point(n))<<std::endl;
            u=1.0;
            l=0.0;
        }
    }
}




template <class Zonotope, class HPolytope, class VT, class PointList, typename NT>
void get_next_zonoball33(Zonotope &Z, std::vector<HPolytope> &HPolySet,
                         HPolytope &HP2, VT Zs_max, VT Zs_min, PointList randPoints,
                         std::vector<NT> &ratios, NT p_value, NT up_lim, int nu, NT alpha){

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension(), count;
    bool too_few;
    bool print = false;

    NT rad2=0.0;
    NT rad1=0.0, rad;
    NT pnorm, ratio;

    VT Zmed(Zs_max.size());
    NT med, u = 1.0, l = 0.0;

    while (true) {
        med = (u + l) * 0.5;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HP2.set_vec(Zmed);
        too_few = false;

        //check_converg001<Point>(HP2, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, false);
        count = 0;
        for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit) {
            if(HP2.is_in(*pit)==-1) count++;
        }
        ratio = NT(count)/NT(1250.0);
        if(ratio>0.09 && ratio<0.12){
            HPolySet.push_back(HP2);
            ratios.push_back(ratio);
            return;
        } else if(ratio>0.12){
            u=med;
        } else {
            l=med;
        }
        std::cout << "ratio = "<<ratio<< std::endl;
        //if(too_few) {
        //    l = med;
        //} else {
        //    u = med;
        //}
    }
}



template <class ZonoHP,class Vpolytope, class HPolytope, class Parameters, typename NT>
void get_sequence_of_vpoly_hpolys33(Vpolytope &VP, HPolytope &HP, std::vector<HPolytope> &HPolySet,
                                  std::vector<NT> &ratios, int Ntot, int nu,
                                  NT &p_value, NT up_lim, NT alpha, Parameters &var, Parameters &var2,
                                  std::vector<NT> &diams_inter) {

    bool print = var.verbose, too_few=false;
    typedef typename Vpolytope::PolytopePoint Point;
    typedef typename Vpolytope::MT MT;
    typedef typename Vpolytope::VT VT;

    int n = var.n, count;

    //MT V = VP.get_mat().transpose();
    //MT AV = HP.get_mat()*V;
    VT Zs_max = VT::Ones(HP.num_of_hyperplanes());

    NT ratio;
    std::list<Point> randPoints;
    Point q(n);
    std::cout<<"sample = "<<Ntot<<" points"<<std::endl;
    std::cout<<"walk_step = "<<var.walk_steps<<std::endl;
    rand_point_generator(VP, q, Ntot, var.walk_steps, randPoints, var);
    var.TotSteps = var.TotSteps + NT(Ntot);
    std::cout<<"points sampled"<<std::endl;
    //check_converg001<Point>(HP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true);
    count = 0;
    HPolytope HP22 = HP;
    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit) {
        if(HP22.is_in(*pit)==-1) count++;
    }
    ratio = NT(count)/NT(Ntot);
    if (ratio >0.09) {
        ratios.push_back(ratio);
        if(print) std::cout<<"last hpoly and ratio = "<<ratio<<std::endl;
        return;
    }
    HPolytope HP2 = HP;
    //HP.print();
    if(print) std::cout<<"not the last hpoly, ratio = "<<ratio<<std::endl;
    get_next_zonoball33(VP, HPolySet, HP2, Zs_max, HP.get_vec(), randPoints, ratios, p_value, up_lim, nu, alpha);
    if(print) std::cout<<"get first hpoly"<<std::endl;

    ZonoHP ZHP2;
    VT Zs_min = HP.get_vec();
    std::pair<Point, NT> pair_diam, pair_diam2;

    while (true) {

        ZHP2 = ZonoHP(VP,HP2);
        q=Point(n);
        randPoints.clear();
        std::cout<<"computing new diameter"<<std::endl;
        //comp_diam_hpoly_vpoly_inter2(ZonoHP &ZHP, MT V, NT AV, VT b, int m, Parameters &var, std::vector<NT> &diams_inter)
        //comp_diam_hpoly_vpoly_inter2<HPolytope>(ZHP2, V, AV, HP2.get_vec(), HP.num_of_hyperplanes(), var2, diams_inter);
        pair_diam = HP2.ComputeInnerBall_from_origin();
        if (HP2.is_all_positive()){
            std::cout<<"{ANN} all bi positives"<<std::endl;
        } else {
            std::cout<<"{ANN} NOT all bi positives!!"<<std::endl;
        }
        diams_inter.push_back(2.0*std::sqrt(NT(n))*pair_diam.second);
        std::cout<<"[annealing] diameter = "<<diams_inter[diams_inter.size()-1]<<std::endl;
        if (pair_diam.second <= 0.0) {
            throw "unbounded";
        }
        pair_diam2 = HP2.ComputeInnerBall();
        if (pair_diam2.second <= 0.0) {
            std::cout<<"radius_from_LP = "<<pair_diam2.second<<std::endl;
        }
        var.diameter = diams_inter[diams_inter.size()-1];
        std::cout<<"sampling N points from ZHP2, walk_length ="<<var.walk_steps<<std::endl;
        rand_point_generator(ZHP2, q, Ntot, var.walk_steps, randPoints, var);
        var.TotSteps = var.TotSteps + NT(Ntot);
        std::cout<<"Checking convergence"<<std::endl;
        check_converg001<Point>(HP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true);
        if (ratio >0.09) {
            ratios.push_back(ratio);
            if(print) std::cout<<"number of hpolys = "<<HPolySet.size()<<std::endl;
            return;
        }
        std::cout<<"getting a new hpoly"<<std::endl;
        get_next_zonoball33(VP, HPolySet, HP2, HP2.get_vec(), Zs_min, randPoints, ratios, p_value, up_lim, nu, alpha);
        if(print) std::cout<<"got hpoly"<<std::endl;
    }


}


#endif
