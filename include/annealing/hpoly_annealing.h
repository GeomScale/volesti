// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef HPOLY_ANNEALING_H
#define HPOLY_ANNEALING_H


template <class Hpolytope, class ZonoHP, class MT, class VT, class Parameters, typename NT>
void comp_diam_hpoly_zono_inter2(ZonoHP &ZHP, MT G, MT AG, VT b, Parameters &var, std::vector<NT> &diams_inter) {

    typedef typename Hpolytope::PolytopePoint Point;
    //typedef typename ZonoHP::NT NT;

    int k = G.cols(), d= ZHP.dimension();

    MT eyes1(k, 2*k);
    eyes1 << MT::Identity(k,k), -1.0*MT::Identity(k,k);

    MT M1(k, 4*k);

    M1 << AG.transpose(), eyes1;
    MT M = M1.transpose();

    VT bb(4*k);
    for (int i = 0; i < 4*k; ++i) bb(i) = (i < 2*k) ? b(i) : 1.0;

    Hpolytope HP;
    HP.init(d, M, bb);

    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    boundary_rand_point_generator(HP, InnerBall.first, 2*d*d, 1, randPoints, var);

    typename std::list<Point>::iterator rpit=randPoints.begin();
    NT max_norm = 0.0, iter_norm;
    for ( ; rpit!=randPoints.end(); rpit++) {
        iter_norm = (G*Eigen::Map<VT>(&(*rpit).get_coeffs()[0], (*rpit).dimension())).norm();
        if (iter_norm > max_norm) max_norm = iter_norm;
    }
    diams_inter.push_back(2.0 * max_norm);

}

template <class VPolytope, class HPolytope, typename NT, class Parameters>
void get_first_poly(VPolytope &VP, HPolytope &HP, NT lb, NT &up_lim, NT &ratio, Parameters &var){

    typedef typename VPolytope::PolytopePoint Point;
    typedef typename VPolytope::MT MT;
    typedef typename VPolytope::VT VT;

    VT Zs_max = HP.get_vec();
    VT Zs_min = VT::Zero(HP.num_of_hyperplanes());
    HPolytope HPiter=HP;

    int n = VP.dimension(), m = Zs_min.size();
    int N = 1200;
    Point q(n);
    bool too_few, print = false;
    std::list<Point> randPoints;

    NT l=0.0, u=1.0, med;
    VT  Zmed(m);
    int count =0;
    Parameters variter = var;
    while(true) {

        count++;
        q=Point(n);
        med = (u + l) * 0.5;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HPiter.set_vec(Zmed);

        randPoints.clear();
        rand_point_generator(HPiter, q, 1200, 10+2*n, randPoints, variter);
        var.MemLps = var.MemLps + 1200.0;
        too_few = false;

        if(check_converg001<Point>(VP, randPoints, lb, up_lim, too_few, ratio, 10, 0.2, true, false)) {
            HP.set_vec(Zmed);
            return;
        }

        if (too_few) {
            u = med;
        } else {
            l = med;
        }
        if(med>0.9) {
            HP.set_vec(Zmed);
            return;
        }
        if(u-l<0.00000000001) {
            u=1.0;
            l=0.0;
        }
    }
}



template <class Polytope, class HPolytope, class VT, typename NT, class Parameters>
void get_hdelta(Polytope &P, HPolytope &HP, VT &Zs_max_gl, NT lb, NT &up_lim, NT &ratio, Parameters &var){

    typedef typename Polytope::PolytopePoint Point;
    typedef typename Polytope::MT MT;
    MT G = P.get_mat().transpose();
    MT A = HP.get_mat();
    int kk = G.cols();
    VT Zs_max = (A*G).cwiseAbs().rowwise().sum();
    Zs_max_gl = Zs_max;
    VT Zs_min = HP.get_vec();
    VT b = HP.get_vec();
    VT b2 = b;
    HPolytope HPiter=HP;

    int n = P.dimension(), m = Zs_max_gl.size();
    int N = 1200;
    Point q(n);
    bool too_few, print = false;
    std::list<Point> randPoints;

    NT l=0.0, u=1.0, med;
    VT  Zmed(m);
    int count =0;
    Parameters variter = var;
    while(true) {

        count++;
        q=Point(n);
        med = (u + l) * 0.5;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HPiter.set_vec(Zmed);

        randPoints.clear();
        rand_point_generator(HPiter, q, 1200, 10+2*n, randPoints, variter);
        too_few = false;

        if(check_converg001<Point>(P, randPoints, lb, up_lim, too_few, ratio, 10, 0.2, true, false)) {
            HP.set_vec(Zmed);
            return;
        }

        if (too_few) {
            u = med;
        } else {
            l = med;
        }
        if(med>0.9) {
            HP.set_vec(Zmed);
            return;
        }
        if(u-l<0.00000000001) {
            u=1.0;
            l=0.0;
        }
    }
}


template <class Zonotope, class HPolytope, class VT, class PointList, typename NT>
void get_next_zonoball22(Zonotope &Z, std::vector<HPolytope> &HPolySet,
                         HPolytope &HP2, VT Zs_max, VT Zs_min, PointList randPoints,
                        std::vector<NT> &ratios, NT p_value, NT up_lim, int nu, NT alpha){

    typedef typename Zonotope::PolytopePoint Point;
    int n = Z.dimension();
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

        if(check_converg001<Point>(HP2, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, false)){
            HPolySet.push_back(HP2);
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

template <class ZonoHP,class Zonotope, class HPolytope, class VT, class Parameters, typename NT>
void get_sequence_of_zonopolys(Zonotope &Z, HPolytope &HP, std::vector<HPolytope> &HPolySet,
                               VT Zs_max, std::vector<NT> &ratios, int Ntot, int nu,
                               NT &p_value, NT up_lim, NT alpha, Parameters &var, Parameters &var2,
                               std::vector<NT> &diams_inter) {

    bool print = var.verbose, too_few=false;
    typedef typename Zonotope::PolytopePoint Point;
    typedef typename Zonotope::MT MT;

    int n = var.n;
    MT G = Z.get_mat().transpose(), AG = HP.get_mat()*G;
    NT ratio;
    std::list<Point> randPoints;
    Point q(n);

    rand_point_generator(Z, q, Ntot, var.walk_steps, randPoints, var);
    HPolytope HP2 = HP;
    if (check_converg001<Point>(HP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true)) {
        ratios.push_back(ratio);
        return;
    }
    get_next_zonoball22(Z, HPolySet, HP2, Zs_max, HP.get_vec(), randPoints, ratios, p_value, up_lim, nu, alpha);

    ZonoHP ZHP2;
    VT Zs_min = HP.get_vec();

    while (true) {

        ZHP2 = ZonoHP(Z,HP2);
        q=Point(n);
        randPoints.clear();
        comp_diam_hpoly_zono_inter2<HPolytope>(ZHP2, G, AG, HP2.get_vec(), var2, diams_inter);
        var.diameter = diams_inter[diams_inter.size()-1];
        rand_point_generator(ZHP2, q, Ntot, var.walk_steps, randPoints, var);
        if (check_converg001<Point>(HP, randPoints, p_value, up_lim, too_few, ratio, nu, alpha, false, true)) {
            ratios.push_back(ratio);
            return;
        }
        get_next_zonoball22(Z, HPolySet, HP2, HP2.get_vec(), Zs_min, randPoints, ratios, p_value, up_lim, nu, alpha);
    }
}


#endif
