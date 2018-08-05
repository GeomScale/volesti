// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
#include <Rcpp.h>
#include <RcppEigen.h>
#include "volume.h"
//#include "../../external/LPsolve/solve_lp.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double vol_R (Rcpp::NumericMatrix A, int walk_len ,double e, Rcpp::NumericVector Chebychev, bool annealing, int win_len,
             int N, double C, double ratio, double frac, bool ball_walk, double delta, bool coord, bool rounding, bool verbose) {

    int n, nexp=1, n_threads=1,i,j;
    NT exactvol(-1.0);
    bool rand_only=false,
	 round_only=false,
	 file=false, 
	 round=rounding, 
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
         birk=false,
         rotate=false,
         experiments=true,
         coordinate=coord;
    HPolytope<NT> P;

    n=A.ncol()-1;
    int m=A.nrow()-1;
    
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);

    if(verbose){
		std::cout<<n<<" "<<m<<std::endl;
		std::cout<<"rnum is: "<<rnum<<std::endl; 
	}
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    
    vars var(rnum,n,walk_len,n_threads,0.0,0.0,0,0.0,0,rng,urdist,urdist1,verbose,rand_only,round,NN,birk,coordinate);
    std::vector<std::vector<NT> > Pin(m+1, std::vector<NT>(n+1));
    std::vector<NT> bin(m);

    for (i=0; i<m+1; i++){
        for(j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
        }
    }
    P.init(Pin);
    //Compute chebychev ball//
    std::pair<Point,NT> CheBall;
    if(Chebychev.size()!=P.dimension()+1){
        CheBall = P.chebyshev_center();
    }else{
        std::vector<NT> temp_p;
        for (int j=0; j<P.dimension(); j++){
          temp_p.push_back(Chebychev[j]);
        }
        Point xc( P.dimension() , temp_p.begin() , temp_p.end() );
        NT radius = Chebychev[P.dimension()];
        CheBall.first = xc; CheBall.second = radius;
    }

    NT vol;
    if (annealing) {
        vars var2(rnum, n, 10 + n / 10, n_threads, 0.0, e, 0, 0.0, 0, rng,
                  urdist, urdist1, verbose, rand_only, round, NN, birk, coordinate);
        vars_g var1(n, walk_len, N, win_len, 1, e, CheBall.second, rng, C, frac, ratio, delta, false, verbose, rand_only, round,
                    NN, birk, ball_walk, coordinate);
        vol = volume_gaussian_annealing(P, var1, var2, CheBall);
        if(verbose) std::cout<<"volume computed = "<<vol<<std::endl;
    } else {
        vol = volume(P, var, var, CheBall);
    }
    
    return ((double)vol);
}
