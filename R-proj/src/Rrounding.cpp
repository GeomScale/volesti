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
#include "rounding.h"
#include "extractMatPoly.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix rounding (Rcpp::NumericMatrix A, unsigned int walk_len, bool coord,
                              bool ball_walk, double delta, bool Vpoly, bool Zono) {

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

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose=false,
            coordinate=coord;

    unsigned int m = A.nrow() - 1;
    unsigned int n = A.ncol() - 1;
    unsigned int rnum = std::pow(1.0,-2.0) * 400 * n * std::log(n);
    std::vector <std::vector<NT> > Pin(m + 1, std::vector<NT>(n + 1));

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
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericVector vec(n + 1);
    if (Zono) {
        InnerBall = ZP.ComputeInnerBall();
    } else if (!Vpoly) {
        InnerBall = HP.ComputeInnerBall();
    } else {
        InnerBall = VP.ComputeInnerBall();
    }


    Rcpp::NumericMatrix Mat;
    if (ball_walk) {
        delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
    }
    // initialization
    vars<NT, RNGType> var(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                          delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
    std::pair <NT, NT> round_res;
    if (Zono) {
        round_res = rounding_min_ellipsoid(ZP, InnerBall, var);
        Mat = extractMatPoly(ZP);
    } else if (!Vpoly) {
        round_res = rounding_min_ellipsoid(HP, InnerBall, var);
        Mat = extractMatPoly(HP);
    } else {
        round_res = rounding_min_ellipsoid(VP, InnerBall, var);
        Mat = extractMatPoly(VP);
    }
    // store rounding value and the ratio between min and max axe in the first row
    // the matrix is in ine format so the first row is useless and is going to be removed by R function modifyMat()
    Mat(0,0) = round_res.first;
    Mat(0,1) = round_res.second;
    return Mat;

}
