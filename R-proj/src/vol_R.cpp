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
double vol_R(Rcpp::NumericMatrix A, int W ,double e, Rcpp::NumericVector Chebychev, bool coord, bool rounding, bool V){
    
    int n, nexp=1, n_threads=1,i,j;
    int walk_len;//to be defined after n
    double exactvol(-1.0);
    bool verbose=true, 
	 rand_only=false, 
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
    Polytope<double> P;
         
    walk_len=W;
    n=A.ncol()-1;
    int m=A.nrow()-1;
    
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    if(!V){
		verbose=false;
	}
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
    std::vector<std::vector<double> > Pin(m+1, std::vector<double>(n+1));
    std::vector<double> bin(m);
    
    for (i=0; i<m+1; i++){
        //bin[i]=b[i];
        for(j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
            if(verbose){
				std::cout<<Pin[i][j]<<" ";
			}
        }
        if(verbose){
			std::cout<<"\n";
		}
    }
    P.init(Pin);
    //Compute chebychev ball//
    std::pair<Point,double> CheBall;
    if(Chebychev.size()!=P.dimension()+1){
        CheBall = solveLP(P);
    }else{
        std::vector<double> temp_p;
        for (int j=0; j<P.dimension(); j++){
          temp_p.push_back(Chebychev[j]);
        }
        Point xc( P.dimension() , temp_p.begin() , temp_p.end() );
        NT radius = Chebychev[P.dimension()];
        CheBall.first = xc; CheBall.second = radius;
    }
    
    Polytope<double> P_to_test(P);
    
    NT vol = volume(P,var,var,CheBall);
    
    
    return vol;
}
