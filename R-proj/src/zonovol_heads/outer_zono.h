// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

//#include <Rcpp.h>

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


template <class Point, class Polytope, class VT, class MT, typename NT, class PointList>
void get_delta(Polytope &P, VT &l, VT &u, MT &sigma, Rcpp::Function mvrandn, MT G, NT &var, NT &delta, NT &up_lim, NT &ratio, PointList &randPoints){

    NT delta1 = 0.0, delta2;
    if (delta!=0.0){
        delta2 = delta;
    } else {
        delta2 = 0.5;
    }

    if(up_lim==0.0){
        up_lim=0.6;
    }
    int n = P.dimension(), m = P.num_of_vertices();
    int N = 1200;
    if (var==0.0) {
        var = 100.0*m;
    }
    //var = 100.0*n;
    VT l2(m), u2(m);
    MT sigma2(m,m);
    MT sample;

    bool done = false, too_few = false;

    while(true) {

        l2 = l - VT::Ones(m) * delta2;
        u2 = u + VT::Ones(m) * delta2;
        sigma2 = var * sigma;

        sample = sampleTr(l2, u2, sigma2, N, mvrandn, G);
        for (int i = 0; i < N; ++i) {
            Point p(n, typename std::vector<NT>::iterator(sample.col(i).data()), typename std::vector<NT>::iterator(sample.col(i).data() + n));
            randPoints.push_back(p);
        }

        check_converg(P, randPoints, 0.1, done, too_few, ratio, up_lim, true);
        //std::cout<<"ratio = "<<ratio<<std::endl;
        //std::cout<<"delta2 = "<<delta2<<std::endl;

        if(done) {
            delta = delta2;
            return;
        }

        if (too_few) {
            break;
        }

        delta1 = delta2;
        delta2 = 2*delta2;
        var = 2*var;
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

        sample = sampleTr(l2, u2, sigma2, N, mvrandn, G);
        for (int i = 0; i < N; ++i) {
            Point p(n, typename std::vector<NT>::iterator(sample.col(i).data()), typename std::vector<NT>::iterator(sample.col(i).data() + n));
            randPoints.push_back(p);
        }

        check_converg(P, randPoints, 0.1, done, too_few, ratio, up_lim, true);
        //std::cout<<"ratio = "<<ratio<<std::endl;
        //std::cout<<"delta_med = "<<delta_med<<std::endl;

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

}