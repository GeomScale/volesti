#include <Rcpp.h>
#include <vol_rand.h>

double volume(NumericMatrix A, NumericVector b, int W, double e){
    
    int n, nexp=1, n_threads=1,i,j;
	int walk_len;//to be defined after n
    double exactvol(-1.0);
    bool verbose=false, 
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
    double vol=0.0;
    stdHPolytope<double> P;
         
    walk_len=W;
    n=A.ncol();
    int m=A.nrow();
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    
    vars var(rnum,n,walk_len,verbose,rand_only,round,NN,birk,coordinate);
    vector<vector<double>> Pin(m, vector<double>(n));
    vector<double> bin(m);
    
    for (i=0; i<m; i++){
        bin[i]=b[i];
        for(j=0; j<n; j++){
            Pin[i][j]=A(i,j);
        }
    }
    P.init(n,Pin,bin);
    
    
    return vol;
}
