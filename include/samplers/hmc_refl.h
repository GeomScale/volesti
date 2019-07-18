// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_REFL_H
#define HMC_REFL_H

template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_gaussian_ref(Polytope &P, Point &p, NT &a, int N, int walk_step, PointList &randPoints, NT radius = -1.0,
                      NT R = -1.0) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::normal_distribution<> rdist(0, 1);

    MT A = P.get_mat(), At = A.transpose();
    unsigned int d = P.dimension(), m = A.rows(), facet;
    VT b = P.get_vec(), x0(d), v0(d);
    NT t, sumt;

    if (R < 0.0) R = 1.0;
    if (radius > 0.0) {

        r = R * radius;
        A = A * (MT::Identity(d, d) * r);

    }
    for (int j = 0; j < d; ++j) x(j) = (1.0 / r) * p[j];

    for (int i = 0; i < N; ++i) {

        T = urdist(rng) * L;
        sumt = 0.0;
        for (int i = 0; i < d; ++i) v0(i) = rdist(rng);
        if (urdist(rng) > 0.5) v0 = -v0;

        for (int l = 0; l < walk_step; ++l) {

            while(sumt<T) {

                compute_boundary_intersection(A, b, x0, v0, t, facet);

                if (t < T) break;

                sumt += t;
                update_position(x0, v0, t);
            }

        }
        for (int k = 0; k < d; ++k) p.set_coord(k, r * x0(k));
        randPoints.push_back(p);

    }



}


#endif
