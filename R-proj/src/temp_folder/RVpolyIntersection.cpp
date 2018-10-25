// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include "volume.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

double RVpolyIntersection(Rcpp::NumericMatrix V1, Rcpp::NumericMatrix V2,
                          unsigned int walk_len, double e, Rcpp::NumericVector InnerBall,
                          bool CG, unsigned int win_len, unsigned int N, double C,
                          double ratio, double frac, bool ball_walk, double delta,
                          bool coord, bool rounding) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef IntersectionOfVpoly<Vpolytope> InterVP;
    //unsigned int n_threads=1,i,j;

    bool rand_only=false,
            NN=false,
            birk=false,
            verbose =false,
            coordinate=coord;
    unsigned int n_threads=1;
    Vpolytope VP1;
    Vpolytope VP2;

    unsigned int m1=V1.nrow() - 1;
    unsigned int n=V1.ncol() - 1;
    unsigned int m2=V2.nrow() - 1;
    unsigned int rnum = std::pow(e,-2) * 400 * n * std::log(n);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::vector<std::vector<NT> > Pin1(m1+1, std::vector<NT>(n+1));
    std::vector<std::vector<NT> > Pin2(m2+1, std::vector<NT>(n+1));

    for (unsigned int i=0; i<m1+1; i++){
        for(unsigned int j=0; j<n+1; j++){
            Pin1[i][j]=V1(i,j);
        }
    }
    for (unsigned int i=0; i<m2+1; i++){
        for(unsigned int j=0; j<n+1; j++){
            Pin2[i][j]=V2(i,j);
        }
    }
    // construct polytopes
    VP1.init(Pin1);
    VP2.init(Pin2);
    InterVP VPcVP(VP1, VP2);

    std::pair<Point,NT> InnerB;

    if(InnerBall.size()==n+1) { //if it is given as an input
        // store internal point hat is given as input
        std::vector<NT> temp_p;
        for (unsigned int j=0; j<n; j++){
            temp_p.push_back(InnerBall[j]);
        }
        InnerB.first = Point( n , temp_p.begin() , temp_p.end() );
        // store the radius of the internal ball that is given as input
        InnerB.second = InnerBall[n];
    } else {
        // no internal ball or point is given as input
        InnerB = VPcVP.ComputeInnerBall();
        if (InnerB.second < 0.0) {
            return -1.0;
        }
    }


    // initialization
    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,0.0,0.0,0,0.0,0, InnerB.second,rng,urdist,urdist1,
                          delta,verbose,rand_only,rounding,NN,birk,ball_walk,coordinate);
    NT vol;
    if (CG) {
        vars<NT, RNGType> var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                               urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
        vars_g<NT, RNGType> var1(n, walk_len, N, win_len, 1, e, InnerB.second, rng, C, frac, ratio, delta, false, verbose,
                                 rand_only, rounding, NN, birk, ball_walk, coordinate);
        vol = volume_gaussian_annealing(VPcVP, var1, var2, InnerB);

    } else {
        vol = volume(VPcVP, var, var, InnerB);
    }

    return vol;

}