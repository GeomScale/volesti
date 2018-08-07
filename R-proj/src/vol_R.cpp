// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
#include <Rcpp.h>
#include <RcppEigen.h>
#include "use_double.h"
#include "volume.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double vol_R (Rcpp::NumericMatrix A, int walk_len, double e, Rcpp::NumericVector Chebychev, bool annealing, int win_len,
             int N, double C, double ratio, double frac, bool ball_walk, double delta, bool Vpoly, bool coord, bool rounding, bool verbose) {

    int n, nexp=1, n_threads=1,i,j;
    NT exactvol(-1.0);
    bool rand_only=false,
	 round_only=false,
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

    n=A.ncol()-1;
    int m=A.nrow()-1;
    
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);

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

    //Compute chebychev ball//
    std::pair<Point,NT> CheBall;
    if(Chebychev.size()!=P.dimension()+1){ //if it is not given as an input
        if (!Vpoly) {
            CheBall = P.chebyshev_center();
        } else {
            CheBall = VP.chebyshev_center();
        }
    }else{ // if it is given as an input
        std::vector<NT> temp_p;
        for (int j=0; j<n; j++){
          temp_p.push_back(Chebychev[j]);
        }
        Point xc( n , temp_p.begin() , temp_p.end() );
        NT radius = Chebychev[n];
        CheBall.first = xc; CheBall.second = radius;
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
    
    return ((double)vol);
}


Rcpp::NumericMatrix sample_R(Rcpp::NumericMatrix A, int walk_len, int rnum,  Rcpp::NumericVector internal_pont, bool gaussian,
                             double variance, bool ball_walk, double delta, bool Vpoly, bool coord, bool verbose){

    bool rand_only=false, NN=false, birk=false;

    HPolytope<NT> P;
    VPolytope<NT> VP;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    //boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    int n=A.cols()-1, m = A.rows()-1;
    NT a;
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
    std::pair CheBall;
    if (!Vpoly) {
        CheBall = P.chebyshev_center();
    } else {
        CheBall = VP.chebyshev_center();
    }
    if (internal_pont.size() != d){// if it is given as an input
        std::vector<NT> temp_p;
        for (int j=0; j<n; j++){
            temp_p.push_back(internal_pont[j]);
        }
        //Point q( n , temp_p.begin() , temp_p.end() );
        CheBall.first = Point( n , temp_p.begin() , temp_p.end() );
    }
    Point p = CheBall.first;

    std::list<Point> randPoints;
    vars var(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,CheBall.second,rng,urdist,urdist1,
             delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
    if (ball_walk) {
        if (delta<0.0){
            if(!gaussian) {
                var.delta = var.delta = 4.0 * radius / std::sqrt(NT(n));
            } else {
                a = 1.0 / (2.0 * variance);
                var.delta = var.delta = 4.0 * radius / std::sqrt(std::max(FT(1.0), a) * FT(n));
            }
        }
    }
    if (!Vpoly) {
        rand_point_generator(P, p, 1, 50 * n, randPoints, var);
    } else {
        rand_point_generator(VP, p, 1, 50 * n, randPoints, var);
    }
    randPoints.clear();
    if (!gaussian){
        if (!Vpoly) {
            rand_point_generator(P, p, rnum, walk_len, randPoints, var);
        } else {
            rand_point_generator(VP, p, rnum, walk_len, randPoints, var);
        }
    } else {
        vars_g var2(n, walk_len, 0, 0, 1, 0, CheBall.second, rng, 0, 0, 0, delta, false, verbose,
                    rand_only, false, NN, birk, ball_walk, coord);
        rand_gaussian_point_generator(P, p, rnum, walk_len, randPoints, a, var2);
    }

    Rcpp::NumericMatrix PointSet(d,rnum);

    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<FT>::iterator qit;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            PointSet(i,j)=*qit;
        }
    }

    return PointSet;

}

