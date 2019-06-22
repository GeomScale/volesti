// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_RK4_H
#define HMC_RK4_H


template <class RNGType, class VT, class MT, typename NT>
void hmc_logbarrier_rk4(Polytope &P, Point &p, PointList &randPoints, NT &a, NT L=0.0, int N) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::normal_distribution<> rdist(0,1);

    MT A = P.get_mat();
    MT At = A.transpose();
    unsigned int d = P.dimension();
    unsigned int m = A.rows();
    VT b = P.get_vec(), s0(m), sv0(m), Y(2*d), Y05(2*d), mi(2*d), ms = VT::Zero(2*d), v0(d), x0(d);
    NT sumh, h, har = 0.1, T;

    for (int i = 0; i < d; ++i) x0(i) = p[i];

    for (int i = 0; i < N; ++i) {

        T = urdist(rng) * L;
        for (int i = 0; i < d; ++i) v0(i) = rdist(rng);
        if (rdist(rng>0.5)) v0 = -v0;

        Y.segment(0, d-1) = x0; Y.segment(d, 2*d-1) = v0;
        sumh = 0.0;

        while (sumh < T) {

            s0 = A * Y.segment(0,d-1);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d-1) = Y.segment(d, 2*d-1);
            mi.segment(d, 2*d-1) = -At * s0;
            ms += mi;

            Y05 = Y + (0.5 * h) * mi;
            s0 = A * Y05.segment(0,d-1);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d-1) = Y.segment(d, 2*d-1);
            mi.segment(d, 2*d-1) = -At * s0;
            ms += 2.0 * mi;

            Y05 = Y + (0.5 * h) * mi;
            s0 = A * Y05.segment(0,d-1);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d-1) = Y.segment(d, 2*d-1);
            mi.segment(d, 2*d-1) = -At * s0;
            ms += 2.0 * mi;

            Y05 = Y + h * mi;
            s0 = A * Y05.segment(0,d-1);
            for (int j = 0; j < m; ++j) s0(j) = a / (b(j) - s0(j));
            mi.segment(0, d-1) = Y.segment(d, 2*d-1);
            mi.segment(d, 2*d-1) = -At * s0;
            ms += mi;

            Y05 = Y + (h/0.6) * ms;
            s0 = A * Y05.segment(0,d-1);
            if ((s0 - b).maxCoeff() > 0.0){
                h = h*0.5;
                continue;
            }
            Y = Y05;
            sumh += h;
            if (T-sumh > har) {
                h = har;
            } else {
                h = T - har;
            }

        }
        for (int k = 0; k < ; ++k) {
            p[k] = Y(k);
        }
        randPoints.push_back(p);

    }

}


#endif
