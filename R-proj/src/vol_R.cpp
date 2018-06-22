#include <Rcpp.h>
#include <RcppEigen.h>
#include "../../include/comp_vol.h"
//#include "../../include/LPsolve/solve_lp.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double vol_R(Rcpp::NumericMatrix A, int W ,double e, Rcpp::NumericVector C){
    
    int n, nexp=1, n_threads=1,i,j;
	int walk_len;//to be defined after n
    double exactvol(-1.0);
    bool verbose=true, 
	 rand_only=false, 
	 round_only=false,
	 file=false, 
	 round=false, 
	 NN=false,
	 user_walk_len=false,
	 linear_extensions=false,
         birk=false,
         rotate=false,
         experiments=true,
         coordinate=true;
    //double vol=0.0;
    stdHPolytope<double> P;
         
    walk_len=W;
    //double e=0.3;
    n=A.ncol()-1;
    int m=A.nrow()-1;
    std::cout<<n<<" "<<m<<std::endl;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    std::cout<<"rnum is: "<<rnum<<std::endl;  
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    std::normal_distribution<double> rdist(0.0,1.0); 
    std::uniform_real_distribution<double> urdist(0.0,1.0);
    std::uniform_real_distribution<double> urdist1(-1.0,1.0);
    
    vars var(rnum,n,walk_len,n_threads,0.0,0.0,0,0.0,0,rng,urdist,urdist1,verbose,rand_only,round,NN,birk,coordinate);
    std::vector<std::vector<double>> Pin(m+1, std::vector<double>(n+1));
    std::vector<double> bin(m);
    
    for (i=0; i<m+1; i++){
        //bin[i]=b[i];
        for(j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
            std::cout<<Pin[i][j]<<" ";
        }
        std::cout<<"\n";
    }
    P.init(Pin);
    //Compute chebychev ball//
	std::pair<Point,double> CheBall;// = solveLP(P.get_matrix(), P.dimension());
    std::vector<double> temp_p;
    for (int j=0; j<P.dimension(); j++){
		temp_p.push_back(C[j]);
	}
	Point xc( P.dimension() , temp_p.begin() , temp_p.end() );
	double radius = C[P.dimension()];
	CheBall.first = xc; CheBall.second = radius;
    
    
    stdHPolytope<double> P_to_test(P);
    
    NT Chebtime;
    
    NT vol = volume1_reuse2(P,var,var,CheBall,Chebtime);
    
    
    return vol;
}
