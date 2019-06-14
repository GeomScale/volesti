// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_HEUR_SAM_H
#define HMC_HEUR_SAM_H


template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_logbarrier(Polytope &P, Point &p, PointList randPoints, NT &a, int N) {

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

    for (int i = 0; i < N; ++i) {

        get_next_hmc_barrier<RNGType>(x0, s0, sv0, A, b, a);

        for (int i = 0; i < d; ++i) {
            v0(i) = rdist(rng);
            p[i] = x0(i);
        }
        randPoints.push_back(p);
        s0 = A*x0;
        sv0 = A*v0;
    }

}


template <class RNGType, class VT, class MT, typename NT>
void get_next_hmc_barrier<RNGType>(VT &x0, VT &s0, VT &sv0, MT &M, VT &b, MT &pinvA, NT &a) {

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::random::uniform_real_distribution<> urdist(0, 1);

    NT tpos = std::numeric_limits<NT>::max(), tminus = std::numeric_limits<NT>::lowest(), t1, t2;

    unsigned int m = s0.size();
    unsigned int d = x0.size();

    VT abcs(m, 3);
    abcs.col(2) = sv0;
    abcs.col(1) = s0;

    VT bt = b - abcs.col(2);
    abcs.col(0) = (0.5 * a) * M * (bt.pow(-1.0));

    abcs.col(2) -= b;
    VT deltas = abcs.col(1).pow(2.0) - 4.0 * abcs.col(0) * abcs.col(2);

    for (int i = 0; i < m; ++i) {

        if (deltas(i) < 0.0) continue;

        t1 = (-abcs(i, 1) - std::sqrt(deltas(i))) / (2.0 * abcs(i, 0));
        if (t1 > 0 && t1 < tpos) tpos = t1;
        if (t1 < 0 && t1 > tminus) tminus = t1;

        t2 = (-abc(i, 1) + std::sqrt(deltas(i))) / (2.0 * abc(i, 0));
        if (t2 > 0 && t2 < tpos) tpos = t2;
        if (t2 < 0 && t2 > tminus) tminus = t2;

    }

    NT trand = urdist(rng) * (tpos - tminus) + tminus;
    VT tvec(3);
    tvec(0) = trand * trand;
    tvec(1) = trand;
    tvec(2) = 1.0;
    x0 = (pinvA * abcs) * tvec;

}


#endif
