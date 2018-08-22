// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
#include <Rcpp.h>
#include <RcppEigen.h>
#include "use_double.h"
#include "volume.h"
#include "extractMatPoly.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix vol_R (Rcpp::NumericMatrix A, int walk_len, double e, Rcpp::NumericVector Chebychev, bool annealing, int win_len,
             int N, double C, double ratio, double frac, bool ball_walk, double delta, bool Vpoly, bool round_only, bool rotate_only, bool sample_only, int numpoints, double variance, bool coord, bool rounding, bool verbose) {


    int nexp=1, n_threads=1,i,j;
    NT exactvol(-1.0);
    bool rand_only=false,
	 file=false,
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
         birk=false,
         rotate=false,
         experiments=true,
         coordinate=coord;
    HPolytope<NT> P;
    VPolytope<NT> VP;

    int m=A.nrow()-1;
    int n=A.ncol()-1;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    Rcpp::NumericMatrix vol_res(1,1);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    //boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::vector<std::vector<NT> > Pin(m+1, std::vector<NT>(n+1));
    std::vector<NT> bin(m);

    for (i=0; i<m+1; i++){
        for(j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
        }
    }
    // construct polytope
    if (!Vpoly) {
        P.init(Pin);
    } else {
        VP.init(Pin);
    }

    if (rotate_only) {
        Rcpp::NumericMatrix Mat;
        if (!Vpoly) {
            rotating(P);
            Mat = extractMatPoly(P);
        } else {
            rotating(VP);
            Mat = extractMatPoly(VP);
        }
        return Mat;
    }

    std::pair<Point,NT> InnerBall;
    //Compute chebychev ball//
    if(Chebychev.size()==n+1 || Chebychev.size()==n) { //if it is given as an input
        if (Chebychev.size()==n) {
            // if only sampling is requested
            // the radius of the inscribed ball is going to be needed for the sampling (radius of ball walk)
            if (!Vpoly) {
                InnerBall = P.ComputeInnerBall();
            } else {
                InnerBall = VP.ComputeInnerBall();
            }
        }
        // store internal point hat is given as input
        std::vector<NT> temp_p;
        for (int j=0; j<n; j++){
            temp_p.push_back(Chebychev[j]);
        }
        InnerBall.first = Point( n , temp_p.begin() , temp_p.end() );
        // store the radius of the internal ball that is given as input
        if (Chebychev.size()==n+1) InnerBall.second = Chebychev[n];
    } else {
        // no internal ball or point is given as input
        if (!Vpoly) {
            InnerBall = P.ComputeInnerBall();
        } else {
            InnerBall = VP.ComputeInnerBall();
        }
    }

    // if only rounding is requested
    if (round_only) {
        Rcpp::NumericMatrix Mat;
        if (ball_walk) {
            delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
        }
        // initialization
        vars var(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                 delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
        std::pair <NT, NT> round_res;
        if (!Vpoly) {
            round_res = rounding_min_ellipsoid(P, InnerBall, var);
            Mat = extractMatPoly(P);
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

    // if only sampling is requested
    if (sample_only){
        std::list<Point> randPoints;
        Point p = InnerBall.first;
        NT a = 1.0 / (2.0 * variance);
        if (ball_walk){
            if(delta<0.0){ // set the radius for the ball walk if is not set by the user
                if(annealing) {
                    delta = 4.0 * InnerBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
                } else {
                    delta = 4.0 * InnerBall.second / std::sqrt(NT(n));
                }
            }
        }
        // initialization
        vars var1(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
                 delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
        vars_g var2(n, walk_len, 0, 0, 1, 0, InnerBall.second, rng, 0, 0, 0, delta, false, verbose,
                    rand_only, false, NN, birk, ball_walk, coord);
        if(!Vpoly) {
            sampling_only(randPoints, P, walk_len, numpoints, annealing, a, p, var1, var2);
        } else {
            sampling_only(randPoints, VP, walk_len, numpoints, annealing, a, p, var1, var2);
        }
        Rcpp::NumericMatrix PointSet(n,numpoints);

        // store the sampled points to the output matrix
        typename std::list<Point>::iterator rpit=randPoints.begin();
        typename std::vector<NT>::iterator qit;
        j = 0;
        for ( ; rpit!=randPoints.end(); rpit++, j++) {
            qit = (*rpit).iter_begin(); i=0;
            for ( ; qit!=(*rpit).iter_end(); qit++, i++){
                PointSet(i,j)=*qit;
            }
        }
        return PointSet;
    }

    // print chebychev ball in verbose mode
    if (verbose) {
        std::cout << "Inner ball center = " << std::endl;
        for (i = 0; i < n; i++) {
            std::cout << InnerBall.first[i] << " ";
        }
        std::cout << "\nradius of inner ball = " << InnerBall.second << std::endl;
    }

    // initialization for volesti
    vars var(rnum,n,walk_len,n_threads,0.0,0.0,0,0.0,0,InnerBall.second,rng,urdist,urdist1,
             delta,verbose,rand_only,rounding,NN,birk,ball_walk,coordinate);
    NT vol;
    if (annealing) {
        // initialization for CV
        vars var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, InnerBall.second, rng,
                  urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
        vars_g var1(n, walk_len, N, win_len, 1, e, InnerBall.second, rng, C, frac, ratio, delta, false, verbose,
                    rand_only, rounding, NN, birk, ball_walk, coordinate);
        if (!Vpoly) { // if the input is a H-polytope
            vol = volume_gaussian_annealing(P, var1, var2, InnerBall);
        } else {  // if the input is a V-polytope
            vol = volume_gaussian_annealing(VP, var1, var2, InnerBall);
        }
        if (verbose) std::cout << "volume computed = " << vol << std::endl;
    } else {
        if (!Vpoly) { // if the input is a H-polytope
            vol = volume(P, var, var, InnerBall);
        } else { // if the input is a V-polytope
            vol = volume(VP, var, var, InnerBall);
        }
    }
    vol_res(0,0)=vol;
    
    return vol_res;
}
