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

    int n=A.ncol()-1;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    //if (sample_only) {
        //return sample_R(A, walk_len, numpoints, Chebychev, annealing, variance, ball_walk, delta, Vpoly, coord, verbose);
    //}
    Rcpp::NumericMatrix vol_res(1,1);

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

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
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

    std::pair<Point,NT> CheBall;
    //Compute chebychev ball//
    if(Chebychev.size()==n+1 || Chebychev.size()==n) { //if it is given as an input
        if (Chebychev.size()==n) {
            if (!Vpoly) {
                CheBall = P.chebyshev_center();
            } else {
                CheBall = VP.chebyshev_center();
            }
        }
        std::vector<NT> temp_p;
        for (int j=0; j<n; j++){
            temp_p.push_back(Chebychev[j]);
        }
        CheBall.first = Point( n , temp_p.begin() , temp_p.end() );
        if (Chebychev.size()==n+1) CheBall.second = Chebychev[n];
    } else {
        if (!Vpoly) {
            CheBall = P.chebyshev_center();
        } else {
            CheBall = VP.chebyshev_center();
        }
    }

    if (round_only) {
        Rcpp::NumericMatrix Mat;
        if (ball_walk) {
            delta = 4.0 * CheBall.second / std::sqrt(NT(n));
        }
        vars var(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,CheBall.second,rng,urdist,urdist1,
                 delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
        if (!Vpoly) {
            rounding_min_ellipsoid(P, CheBall, var);
            Mat = extractMatPoly(VP);
        } else {
            rounding_min_ellipsoid(VP, CheBall, var);
            Mat = extractMatPoly(VP);
        }
        return Mat;
    }

    if (sample_only){
        std::list<Point> randPoints;
        Point p = CheBall.first;
        NT a = 1.0 / (2.0 * variance);
        if (ball_walk){
            if(delta<0.0){
                if(annealing) {
                    delta = 4.0 * CheBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
                } else {
                    delta = 4.0 * CheBall.second / std::sqrt(NT(n));
                }
            }
        }
        vars var1(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,CheBall.second,rng,urdist,urdist1,
                 delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
        vars_g var2(n, walk_len, 0, 0, 1, 0, CheBall.second, rng, 0, 0, 0, delta, false, verbose,
                    rand_only, false, NN, birk, ball_walk, coord);
        if(!Vpoly) {
            sampling_only(randPoints, P, walk_len, numpoints, annealing, a, p, var1, var2);
        } else {
            sampling_only(randPoints, VP, walk_len, numpoints, annealing, a, p, var1, var2);
        }
        Rcpp::NumericMatrix PointSet(n,numpoints);

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
        std::cout << "Chebychev center = " << std::endl;
        for (i = 0; i < n; i++) {
            std::cout << CheBall.first[i] << " ";
        }
        std::cout << "\nradius of chebychev ball = " << CheBall.second << std::endl;
    }

    // initialization
    vars var(rnum,n,walk_len,n_threads,0.0,0.0,0,0.0,0,CheBall.second,rng,urdist,urdist1,
             delta,verbose,rand_only,rounding,NN,birk,ball_walk,coordinate);
    NT vol;
    if (annealing) {
        vars var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, CheBall.second, rng,
                  urdist, urdist1, delta, verbose, rand_only, rounding, NN, birk, ball_walk, coordinate);
        vars_g var1(n, walk_len, N, win_len, 1, e, CheBall.second, rng, C, frac, ratio, delta, false, verbose,
                    rand_only, rounding, NN, birk, ball_walk, coordinate);
        if (!Vpoly) { // if the input is a H-polytope
            vol = volume_gaussian_annealing(P, var1, var2, CheBall);
        } else {  // if the input is a V-polytope
            vol = volume_gaussian_annealing(VP, var1, var2, CheBall);
        }
        if (verbose) std::cout << "volume computed = " << vol << std::endl;
    } else {
        if (!Vpoly) { // if the input is a H-polytope
            vol = volume(P, var, var, CheBall);
        } else { // if the input is a V-polytope
            vol = volume(VP, var, var, CheBall);
        }
    }
    vol_res(0,0)=vol;
    
    return vol_res;
}
