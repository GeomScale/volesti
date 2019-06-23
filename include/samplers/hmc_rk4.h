// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_RK4_H
#define HMC_RK4_H


template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier_rk4(Polytope &P, Point &p, PointList &randPoints, NT &a, int N,  NT radius = -1.0, NT R = -1.0) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    std::cout<<"hello rk4!"<<std::endl;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::normal_distribution<> rdist(0,1);

    MT A = P.get_mat();
    MT At = A.transpose();
    unsigned int d = P.dimension();
    unsigned int m = A.rows();
    VT b = P.get_vec(), s0(m), sv0(m), Y(2*d+1), Y05(2*d+1), mi(2*d+1), ms = VT::Zero(2*d+1), v0(d), x0(d);
    NT sumh, h, har = 0.1, T, r = 1.0, L = 1.0;
    bool check;

    if (R < 0.0) R = 1.0;
    if (radius > 0.0) {
        std::pair<Point, NT> InnerBall = P.ComputeInnerBall();
        r = R * InnerBall.second;
        std::cout<<MT::Identity(d,d)<<"\n"<<r<<std::endl;
        A = A * (MT::Identity(d,d) * r);
        At = A.transpose();
    }

    for (int i = 0; i < d; ++i) x0(i) = (1 / r) * p[i];
    std::cout<<"x0 = "<<x0<<"\n"<<std::endl;
    std::cout<<"A*x0-b = "<<A*x0 - b<<"\n"<<std::endl;

    //for (int l = 0; l < N; ++l) {
   //     randPoints.push_back(p);
//    }
    //return;

    for (int i = 0; i < N; ++i) {

        T = urdist(rng) * L;
        for (int i = 0; i < d; ++i) v0(i) = rdist(rng);
        if (urdist(rng)>0.5) v0 = -v0;

        //std::cout<<"T = "<<T<<std::endl;
        //std::cout<<"v0 = "<<v0<<"\n"<<std::endl;
        Y.segment(0, d) = x0; Y.segment(d, 2*d) = v0;
        //std::cout<<Y.segment(0, d)<<"\n"<<std::endl;
        //std::cout<<Y.segment(d, 2*d)<<"\n"<<std::endl;
        //std::cout<<Y<<"\n"<<std::endl;
        //std::cout<< A * Y.segment(0,d)<<"\n"<<std::endl;
        //std::cout<< A * x0<<"\n"<<std::endl;
        //mi.segment(0, d) = Y.segment(d, 2*d);
        //std::cout<< mi<<"\n"<<std::endl;
        //s0 = A * Y.segment(0,d);
        //mi.segment(d, 2*d) = -At * s0;
        //std::cout<< -At * s0<<"\n"<<std::endl;
        //std::cout<< mi<<"\n"<<std::endl;


        sumh = 0.0;
        if (T > har) {
            h = har;
        } else {
            h = har - T;
        }
        //h = har;

        while (sumh < T) {

            ms = VT::Zero(2*d);
            s0 = A * Y.segment(0,d);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d) = Y.segment(d, 2*d);
            mi.segment(d, 2*d) = -At * s0;
            ms += mi;

            Y05 = Y + (0.5 * h) * mi;
            s0 = A * Y05.segment(0,d);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d) = Y05.segment(d, 2*d);
            mi.segment(d, 2*d) = -At * s0;
            ms += 2.0 * mi;

            Y05 = Y + (0.5 * h) * mi;
            s0 = A * Y05.segment(0,d);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d) = Y05.segment(d, 2*d);
            mi.segment(d, 2*d) = -At * s0;
            ms += 2.0 * mi;

            Y05 = Y + h * mi;
            s0 = A * Y05.segment(0,d);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d) = Y05.segment(d, 2*d);
            mi.segment(d, 2*d) = -At * s0;
            ms += mi;

            Y05 = Y + (h/6.0) * ms;
            s0 = A * Y05.segment(0,d);
            //std::cout<<s0-b<<"\n"<<std::endl;
            std::cout<<" T = "<<T<<", i = "<<i<<", sumh = "<<sumh<<", h = "<<h<<", maxCoeff = "<<(s0 - b).maxCoeff()<<std::endl;
            if ((s0 - b).maxCoeff() > 0.0){
                h = h*0.5;
                //std::cout<<"i = "<<i<<" sumh = "<<sumh<<" h = "<<h<<std::endl;
                continue;
            }
            Y = Y05;
            sumh += h;
            if (T-sumh > har) {
                h = har;
            } else {
                h = T - sumh;
            }
            check = sumh < T;
            std::cout<<"h = "<<h<<" sumh<T = "<<check<<std::endl;

        }
        for (int k = 0; k < d; ++k) {
            p.set_coord(k, r * Y(k));
        }
        std::cout<<"is in P = "<<P.is_in(p)<<std::endl;
        randPoints.push_back(p);

    }

    std::cout<<"hmc_rk4 ended!"<<std::endl;

}


#endif
