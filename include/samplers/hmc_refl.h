// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HMC_REFL_H
#define HMC_REFL_H


template <class VT, class VT2>
void compute_reflection(VT &v0, VT2 nH) {

    nH = nH / nH.norm();
    v0 = v0 - (2.0 * v0.dot(nH)) * nH;

}

template <typename NT>
void compute_tranj_params(NT xi, NT vi, NT &a, NT &omega, NT &C, NT &Phi) {

    omega = std::sqrt(2.0 * a);
    C = std::sqrt(xi*xi + (vi*vi)/(omega * omega));
    Phi = std::atan((-vi) / (xi * omega));
    if (vi < 0.0 && Phi < 0.0) {
        Phi += M_PI;
    } else if (vi > 0.0 && Phi > 0.0) {
        Phi -= M_PI;
    }

}

template <class VT, typename NT>
void update_position_momenta(VT &x0, VT &v0, NT &t, NT &a) {

    NT omega, C, Phi;

    for (int i = 0; i < x0.size(); ++i) {

        compute_tranj_params(x0(i), v0(i) , a, omega, C, Phi);
        x0(i) = C * std::cos(omega * t + Phi);
        v0(i) = -omega * C * std::sin(omega * t + Phi);

    }

}

template <class MT, class VT, typename NT>
bool compute_inter(MT &A, VT &b, VT &x0, VT &v0, NT &a, NT &t, int &facet) {

    bool intersection = false;
    unsigned int m = A.rows();
    NT omega, C, Phi, t1, t2, tmin;
    t = std::numeric_limits<NT>::max();

    for (int i = 0; i < m; ++i) {
        compute_tranj_params(x0.dot(A.row(i)), v0.dot(A.row(i)), a, omega, C, Phi);
        if (C > b(i)) {

            intersection = true;

            t1 = (std::acos(b(i) / C) - Phi) / omega;
            t1 += (t1 < 0.0) ? (2.0*M_PI) / omega : 0.0;

            t2 = (-std::acos(b(i) / C) - Phi) / omega;
            t2 += (t2 < 0.0) ? (2.0*M_PI) / omega : 0.0;

            tmin = std::min(t1, t2);
            if (tmin < t) {
                facet = i;
                t = tmin;
            }

        }
    }

    return intersection;

}


template <class RNGType, class Polytope, class Point, class PointList, typename NT>
void hmc_gaussian_ref(Polytope &P, Point &p, NT &a, int N, int walk_step, PointList &randPoints, NT radius = -1.0,
                      NT R = -1.0) {

    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::normal_distribution<> rdist(0, 1);

    MT A = P.get_mat();
    unsigned int d = P.dimension(), m = A.rows();
    int facet;
    VT b = P.get_vec(), x0(d), v0(d);
    NT t, sumt, r = 1.0, L = 1.0, T;

    if (R < 0.0) R = 1.0;
    if (radius > 0.0) {

        r = R * radius;
        A = A * (MT::Identity(d, d) * r);

    }
    for (int j = 0; j < d; ++j) x0(j) = (1.0 / r) * p[j];

    for (int i = 0; i < N; ++i) {

        for (int l = 0; l < walk_step; ++l) {

            T = urdist(rng) * L;
            sumt = 0.0;
            for (int i = 0; i < d; ++i) v0(i) = rdist(rng);
            if (urdist(rng) > 0.5) v0 = -v0;

            while(T > sumt) {

                if (!compute_inter(A, b, x0, v0, a, t, facet)) {
                    t = T - sumt;
                    update_position_momenta(x0, v0, t, a);
                    break;
                } else {
                    t = (t > T - sumt) ? T - sumt : 0.999 * t;
                    sumt += t;
                    update_position_momenta(x0, v0, t, a);
                    compute_reflection(v0, A.row(facet).transpose());
                }

            }

        }
        for (int k = 0; k < d; ++k) p.set_coord(k, r * x0(k));
        randPoints.push_back(p);

    }



}


#endif
