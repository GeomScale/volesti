// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MIN_ELLIPSOID_ROUNDING_HPP
#define MIN_ELLIPSOID_ROUNDING_HPP

#include <Eigen/Eigen>
#include "khach.h"
#include "sampling/random_point_generators.hpp"
#include "volume/sampling_policies.hpp"

template
<
    typename WalkTypePolicy,
    typename MT,
    typename VT,
    typename Polytope,
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
std::tuple<MT, VT, NT> min_sampling_covering_ellipsoid_rounding(Polytope &P,
                                                                std::pair<Point,NT> &InnerBall,
                                                                const unsigned int &walk_length,
                                                                RandomNumberGenerator &rng)
{
    typedef typename WalkTypePolicy::template Walk
            <
                    Polytope,
                    RandomNumberGenerator
            > WalkType;
    typedef RandomPointGenerator<WalkType> RandomPointGenerator;

    unsigned int d = P.dimension();
    std::list<Point> randPoints; //ds for storing rand points
    NT ratio = 10, round_val = 1.0;
    unsigned int iter = 0, j, i;
    const unsigned int num_of_samples = 10*d;//this is the number of sample points will used to compute min_ellipoid
    PushBackWalkPolicy push_back_policy;

    MT T = MT::Identity(d,d);
    VT shift = VT::Zero(d);

    while (ratio > 6.0 && iter < 3)
    {
        randPoints.clear();
        if (!P.get_points_for_rounding(randPoints))
        {  // If P is a V-polytope then it will store its vertices in randPoints
            // If P is not a V-Polytope or number_of_vertices>20*domension
            // 2. Generate the first random point in P
            // Perform random walk on random point in the Chebychev ball
            Point c = InnerBall.first;
            NT radius = InnerBall.second;
            Point p = GetPointInDsphere<Point>::apply(d, radius, rng);
            p += c;
            RandomPointGenerator::apply(P, p, num_of_samples, walk_length,
                                        randPoints, push_back_policy, rng);
        }

        // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
        MT Ap(d, randPoints.size());
        typename std::list<Point>::iterator rpit=randPoints.begin();

        j = 0;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            for (i=0 ; i<rpit->dimension(); i++){
                Ap(i,j)=double((*rpit)[i]);
            }
        }
        MT Q(d,d); //TODO: remove dependence on ublas and copy to eigen
        VT c2(d);
        size_t w=1000;
        KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

        MT E(d, d);
        VT e(d);

        //Get ellipsoid matrix and center as Eigen objects
        for(i=0; i<d; i++){
            e(i)=NT(c2(i));
            for (j=0; j<d; j++){
                E(i,j)=NT(Q(i,j));
            }
        }

        //Find the smallest and the largest axes of the elliposoid
        Eigen::EigenSolver<MT> eigensolver(E);
        NT rel = std::real(eigensolver.eigenvalues()[0]);
        NT Rel = std::real(eigensolver.eigenvalues()[0]);
        for(i=1; i<d; i++) {
            if (std::real(eigensolver.eigenvalues()[i]) < rel) rel = std::real(eigensolver.eigenvalues()[i]);
            if (std::real(eigensolver.eigenvalues()[i]) > Rel) Rel = std::real(eigensolver.eigenvalues()[i]);
        }

        Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
        MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

        //Shift polytope in order to contain the origin (center of the ellipsoid)
        P.shift(e);

        MT L_1 = L.inverse();
        shift = shift + T * e;
        T = T * L_1.transpose();

        P.linear_transformIt(L_1.transpose());
        InnerBall = P.ComputeInnerBall();
        round_val *= L_1.determinant();
        ratio = Rel / rel;
        iter++;
    }

    std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(round_val));
    return result;
}


#endif // MIN_ELLIPSOID_ROUNDING_HPP
