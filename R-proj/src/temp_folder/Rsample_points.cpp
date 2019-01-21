// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "polytopes.h"
#include "samplers.h"
#include "gaussian_samplers.h"
#include "sample_only.h"
#include "simplex_samplers.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix Rsample_points (Rcpp::NumericMatrix A, unsigned int walk_len, Rcpp::NumericVector InnerPoint,
                                    bool gaussian, bool ball_walk, double delta, bool coord, bool Vpoly, bool Zono, bool sam_simplex,
                                    bool sam_can_simplex, bool sam_arb_simplex, bool sam_ball, bool sam_sphere,
                                    unsigned int numpoints, int dim, double variance) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> zonotope;

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;

    std::list<Point> randPoints;
    Rcpp::NumericMatrix PointSet(dim,numpoints);

    if (sam_ball || sam_sphere) {

        for (unsigned int k = 0; k < numpoints; ++k) {
            if (sam_ball) {
                randPoints.push_back(get_point_in_Dsphere<RNGType , Point >(dim, delta));
            } else {
                randPoints.push_back(get_point_on_Dsphere<RNGType , Point >(dim, delta));
            }
        }

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        unsigned int j = 0, i;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;

    }

    if (sam_simplex || sam_can_simplex) {

        if (sam_simplex) {
            Sam_Unit<NT, RNGType >(dim, numpoints, randPoints);
        } else {
            Sam_Canon_Unit<NT, RNGType >(dim, numpoints, randPoints);
        }

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        unsigned int j = 0, i;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    if (sam_arb_simplex) {
        unsigned int n=A.ncol()-1;
        std::vector<NT> temp_p(n, 0.0);
        typename std::vector<NT>::iterator temp_it;
        std::vector<Point> vec_point;

        for (int k = 1; k < A.nrow(); ++k) {
            temp_it = temp_p.begin();
            for (int l = 1; l < A.ncol(); ++l, ++temp_it) {
                *temp_it = A(k,l);
            }
            vec_point.push_back(Point(n, temp_p.begin(), temp_p.end()));
        }

        Sam_arb_simplex<NT, RNGType>(vec_point.begin(), vec_point.end(), numpoints, randPoints);

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        unsigned int j = 0, i;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose=false,
            coordinate=coord;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    unsigned int m=A.nrow()-1;
    unsigned int n=A.ncol()-1;
    std::vector <std::vector<NT> > Pin(m + 1, std::vector<NT>(n + 1));
    Rcpp::NumericMatrix PointSet2(n,numpoints);

    for (unsigned int i = 0; i < m + 1; i++) {
        for (unsigned int j = 0; j < n + 1; j++) {
            Pin[i][j] = A(i, j);
        }
    }
    // construct polytope
    if (Zono) {
        ZP.init(Pin);
    } else if (!Vpoly) {
        HP.init(Pin);
    } else {
        VP.init(Pin);
    }

    std::pair<Point,NT> InnerBall;
    if (Zono) {
        InnerBall = ZP.ComputeInnerBall();
    } else if (!Vpoly) {
        InnerBall = HP.ComputeInnerBall();
    } else {
        InnerBall = VP.ComputeInnerBall();
    }

    if (InnerPoint.size()==n) {
        std::vector<NT> temp_p;
        for (unsigned int j=0; j<n; j++){
            temp_p.push_back(InnerPoint[j]);
        }
        InnerBall.first = Point( n , temp_p.begin() , temp_p.end() );
    }

    Point p = InnerBall.first;
    NT a = 1.0 / (2.0 * variance);
    if (ball_walk) {
        if (delta < 0.0) { // set the radius for the ball walk if is not set by the user
            if (gaussian) {
                delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
            } else {
                delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
            }
        }
    }
    unsigned int rnum = std::pow(1.0,-2) * 400 * n * std::log(n);
    // initialization
    vars<NT, RNGType> var1(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                           delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
    vars_g<NT, RNGType> var2(n, walk_len, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, false, verbose,
                             rand_only, false, NN, birk, ball_walk, coord);
    if (Zono) {
        sampling_only<Point>(randPoints, ZP, walk_len, numpoints, gaussian, a, p, var1, var2);
    } else if (!Vpoly) {
        sampling_only<Point>(randPoints, HP, walk_len, numpoints, gaussian, a, p, var1, var2);
    } else {
        sampling_only<Point>(randPoints, VP, walk_len, numpoints, gaussian, a, p, var1, var2);
    }

    // store the sampled points to the output matrix
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    unsigned int j = 0, i;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            PointSet2(i,j)=*qit;
        }
    }
    return PointSet2;

}
