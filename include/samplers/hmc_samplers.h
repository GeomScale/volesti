// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_SAMPLERS_H
#define HMC_SAMPLERS_H

#include <cmath>

template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier(Polytope &P, Point &p, PointList randPoints, NT &a, int n, int N) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);

    MT A = P.get_mat();
    VT b = P.get_vec();
    unsigned int d = P.dimension();
    unsigned int m = A.nrows();

    VT v0(d), x0(d);
    for (int i = 0; i < d; ++i) {
        v0(i) = rdist(rng);
        x0(i) = p[i];
    }

    s0 = A*x0;
    sv0 = A*v0;
    MT M = - A*A.transpose();

    VT c(3);
    cj(0) = 0.0; cj(1) = 0.587785; cj(2) = 0.951056; //Chebyshev nodes

    MT AA = MT::Ones(n+1, n+1);
    MT T = MT::Zero(n+1,n+1);
    VT S = VT::Zero(n+1);

    for (int j = 1; j < n+1; ++j) {
        AA.col(j) = cj.pow(j);
        S(j) = NT(j)*cj,pow(j-1);
        if (j>1) T.col(j) = NT(j)*NT(j-1)*cj.pow(j-2);
    }
    AAinv = A.inverse();
    T = T*AAinv;
    S = S*AAinv;

    for (int i = 0; i < N; ++i) {

        get_next_hmc_logbarrier(x0, T, S, M, s0, sv0);

        for (int i = 0; i < d; ++i) {
            v0(i) = rdist(rng);
            p[i] = x0(i);
        }
        randPoints.push_back(p);
        s0 = A*x0;
        sv0 = A*v0;
    }

}




#endif
