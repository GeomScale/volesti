// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

//#include <Rcpp.h>

#ifndef OUTER_ZONO_H
#define OUTER_ZONO_H


template <class Point, class ball, class PointList, typename NT>
void check_converg00001(ball &P, PointList &randPoints, NT p_test, bool &done, bool &too_few,
                      NT &ratio, NT up_lim, bool print) {

    //typedef typename ball::BallPoint Point;
    std::vector<NT> ratios;
    NT countsIn = 0.0;
    int m = randPoints.size()/10;
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
            if (rm+zp*rs<0.09) {
                too_few = true;
                return;
            } else if (rm-zp*rs>0.17) {
                return;
            }
        }

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
            //std::cout<<"ni = "<<ni<<std::endl;
            //std::cout<<"ratio done = "<<ratio<<" var = "<<p_varval<<" "<<(up_lim) + t_value*(ni-1)*(p_varval/std::sqrt(NT(ni)))<<std::endl;
        }
        ratio = p_mval; // for test only
    } else {
        too_few = true;
        ratio = p_mval; // for test only
    }

}



template <class Polytope, class PointList, typename NT>
void check_converg(Polytope &P, PointList &randPoints, NT p_test, bool &done, bool &too_few, NT &ratio, NT up_lim, bool print) {

    typedef typename Polytope::PolytopePoint Point;
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

/*
template <class Point, class Polytope, class VT, class MT, typename NT, class PointList>
void get_delta(Polytope &P, VT &l, VT &u, MT &sigma, Rcpp::Function rtmvnorm, Rcpp::Function mvrandn,
               Rcpp::Function mvNcdf, MT G, NT &var, NT &delta, NT &up_lim, NT &ratio, int Wst, PointList &randPoints){

    NT delta1 = 0.0, delta2;
    if (delta!=0.0){
        delta2 = delta;
    } else {
        delta2 = 0.5;
    }

    if(up_lim==0.0){
        up_lim=0.2;
    }
    int n = P.dimension(), m = P.num_of_vertices();
    int N = 1200;
    if (var==0.0) {
        var = 1000000.0;
    }
    //var = 100.0*n;
    VT l2(m), u2(m);
    MT sigma2(m,m);
    MT sample;
    NT prob;
    sigma2 = var * sigma;
    int kk = G.cols();

    l2 = l - VT::Ones(m) * delta2;
    u2 = u + VT::Ones(m) * delta2;
    prob = test_botev<NT>(l2, u2, sigma2, 10000, mvNcdf);
    //prob=0.1;
    std::cout<<"prob = "<<prob<<std::endl;
    //int ww = kk*kk/10;

    bool done = false, too_few = false;

    while(true) {

        l2 = l - VT::Ones(m) * delta2;
        u2 = u + VT::Ones(m) * delta2;
        sigma2 = var * sigma;

        //prob = test_botev(l2, u2, sigma2, 10000, mvNcdf);


        if(prob>0.001) {
            sample = sampleTr(l2, u2, sigma2, N, mvrandn, G);
        } else {
            sample = sampleTr_gibbs(l2, u2, sigma2, N, Wst, rtmvnorm, G);
        }

        for (int i = 0; i < N; ++i) {
            Point p(n, typename std::vector<NT>::iterator(sample.col(i).data()), typename std::vector<NT>::iterator(sample.col(i).data() + n));
            randPoints.push_back(p);
        }

        done = false;
        too_few = false;
        check_converg(P, randPoints, 0.1, done, too_few, ratio, up_lim, true);
        std::cout<<"ratio = "<<ratio<<std::endl;
        std::cout<<"delta2 = "<<delta2<<std::endl;

        if(done) {
            delta = delta2;
            return;
        }

        if (too_few) {
            break;
        }

        delta1 = delta2;
        delta2 = 2*delta2;
        //var = 2*var;
        randPoints.clear();

    }

    // bisection between delta1 - delta_2
    NT delta_med;
    randPoints.clear();
    while(true) {

        delta_med = (delta1 + delta2) * 0.5;

        l2 = l - VT::Ones(m) * delta_med;
        u2 = u + VT::Ones(m) * delta_med;
        //sigma2 = var * sigma;

        if(prob>0.001) {
            sample = sampleTr(l2, u2, sigma2, N, mvrandn, G);
        } else {
            sample = sampleTr_gibbs(l2, u2, sigma2, N, Wst, rtmvnorm, G);
        }
        //sample = sampleTr(l2, u2, sigma2, N, mvrandn, G);
        for (int i = 0; i < N; ++i) {
            Point p(n, typename std::vector<NT>::iterator(sample.col(i).data()), typename std::vector<NT>::iterator(sample.col(i).data() + n));
            randPoints.push_back(p);
        }

        done = false;
        too_few = false;
        check_converg(P, randPoints, 0.1, done, too_few, ratio, up_lim, true);
        std::cout<<"ratio = "<<ratio<<std::endl;
        std::cout<<"delta_med = "<<delta_med<<std::endl;

        if(done) {
            delta = delta_med;
            return;
        }

        if (too_few) {
            delta2 = delta_med;
        } else {
            delta1 = delta_med;
        }
        randPoints.clear();

    }

}*/

template <class Point, class Polytope, class HPolytope, class VT, typename NT, class PointList, class Parameters>
void get_hdelta(Polytope &P, HPolytope &HP, VT &Zs_max_gl, NT &up_lim, NT &ratio,
                PointList &randPoints, Parameters &var, NT &steps){

    NT delta1 = 0.0;
    NT delta2 = 0.5;
    //typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    MT G = P.get_mat().transpose();
    MT A = HP.get_mat();
    int kk = G.cols();
    VT Zs_max = (A*G).cwiseAbs().rowwise().sum();
    Zs_max_gl = Zs_max;
    VT Zs_min = HP.get_vec();

    //get_maxZ0<NT>(A, G, Zs_max);
    //std::cout<<Zs_max<<"\n"<<std::endl;
    VT b = HP.get_vec();
    VT b2 = b;
    HPolytope HPiter=HP;

    int n = P.dimension(), m = P.num_of_generators();
    int N = 1200;
    if(up_lim==0.0){
        up_lim=0.15;
    }
    Point q(n);
    bool done, too_few, print = var.verbose;
    //std::list<Point> randPoints;
    randPoints.clear();
    steps = 0.0;

    NT l=0.0, u=1.0, med;
    VT  Zmed(2*m);
    randPoints.clear();
    int count =0;
    while(true) {

        count++;
        q=Point(n);
        med = (u + l) * 0.5;
        Zmed = Zs_min + (Zs_max-Zs_min)*med;
        HPiter.set_vec(Zmed);

        rand_point_generator(HPiter, q, 1200, 10+2*n, randPoints, var);
        steps += 1200.0;

        done = false;
        too_few = false;
        check_converg00001<Point>(P, randPoints, 0.1, done, too_few, ratio, up_lim, false);
        if(print) std::cout<<"ratio = "<<ratio<<std::endl;
        if(print) std::cout<<"Z_med = "<<med<<std::endl;

        if(done) {
            std::cout<<"done first Hpoly"<<std::endl;
            //delta = delta2;
            HP.set_vec(Zmed);
            return;
        }

        if (too_few) {
            u = med;
            //break;
        } else {
            l = med;
        }

        //delta1 = delta2;
        //delta2 = 2*delta2;
        //var = 2*var;
        randPoints.clear();
        if(count>80 || med>0.9) {
            HP.set_vec(Zmed);
            return;
        }

    }

    //std::cout<<HPiter.get_mat()<<"\n"<<HPiter.get_vec()<<std::endl;
    /*while(true) {

        q=Point(n);
        b2 = b + VT::Ones(2*m) * delta2;
        HPiter.set_vec(b2);
        //std::cout<<HPiter.get_vec()<<std::endl;
        //randPoints.clear();
        rand_point_generator(HPiter, q, 1200, 10+n/10, randPoints, var);
        steps += 1200.0;

        done = false;
        too_few = false;
        check_converg(P, randPoints, 0.1, done, too_few, ratio, up_lim, true);
        if(print) std::cout<<"ratio = "<<ratio<<std::endl;
        if(print) std::cout<<"delta2 = "<<delta2<<std::endl;

        if(done) {
            delta = delta2;
            HP.set_vec(b2);
            return;
        }

        if (too_few) {
            break;
        }

        delta1 = delta2;
        delta2 = 2*delta2;
        //var = 2*var;
        randPoints.clear();

    }

    // bisection between delta1 - delta_2
    NT delta_med;
    randPoints.clear();
    while(true) {

        delta_med = (delta1 + delta2) * 0.5;

        q=Point(n);
        b2 = b + VT::Ones(2*m) * delta_med;
        HPiter.set_vec(b2);
        //randPoints.clear();
        rand_point_generator(HPiter, q, 1200, 10+n/10, randPoints, var);
        steps += 1200.0;

        done = false;
        too_few = false;
        check_converg(P, randPoints, 0.1, done, too_few, ratio, up_lim, true);
        if(print) std::cout<<"ratio = "<<ratio<<std::endl;
        if(print) std::cout<<"delta_med = "<<delta_med<<std::endl;

        if(done) {
            delta = delta_med;
            HP.set_vec(b2);
            return;
        }

        if (too_few) {
            delta2 = delta_med;
        } else {
            delta1 = delta_med;
        }
        randPoints.clear();

    }*/

}

#endif
