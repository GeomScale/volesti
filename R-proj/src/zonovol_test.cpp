

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppEigen.h>
// [[Rcpp::depends(BH)]]
#include <iterator>
//#include <fstream>
#include <vector>
#include <list>
//#include <algorithm>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "polytopes.h"
//#include "ellipsoids.h"
//#include "ballintersectconvex.h"
//#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
//#include "gaussian_samplers.h"
//#include "gaussian_annealing.h"
#include "zonovol_heads/sampleTruncated.h"
#include "zonovol_heads/annealing_zono.h"
#include "zonovol_heads/cg_zonovol.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double vol_zono (Rcpp::Reference P, double e, Rcpp::Function mvrandn, bool verbose) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    //unsigned int n_threads=1,i,j;

    bool rand_only = false,
            NN = false,
            birk = false;
            //verbose = false;
    unsigned int n_threads = 1;
    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;

    Rcpp::NumericMatrix A = P.field("G");

    unsigned int m=A.nrow()-1;
    unsigned int n=A.ncol()-1;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    std::vector<std::vector<NT> > Pin(m+1, std::vector<NT>(n+1));

    for (unsigned int i=0; i<m+1; i++){
        for(unsigned int j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
        }
    }
    ZP.init(Pin);

    MT G = (Rcpp::as<MT>(A)).transpose();
    MT ps = G.completeOrthogonalDecomposition().pseudoInverse();
    MT sigma = ps*ps.transpose();

    VT l = VT::Zero(m) - VT::Ones(m);
    VT u = VT::Ones(m);

    MT test = sampleTr(l, u, sigma, 100, mvrandn, G);
    std::pair<Point,NT> InnerB = ZP.ComputeInnerBall();
    NT C = 2.0, frac=0.1, ratio = 1.0-1.0/(NT(n));;
    int W = 4*n*n+500, N = 500 * ((int) C) + ((int) (n * n / 2));

    vars<NT, RNGType> var2(1, n, 10 + n / 10, 1, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                           urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                           false);
    vars_g<NT, RNGType> var1(n, 1, N, W, 1, e, InnerB.second, rng, C, frac,
                             ratio, -1.0, false, verbose, false, false, NN, birk,
                             false, false);

    NT vol = cg_volume_zono(ZP, var1, var2, InnerB, mvrandn, sigma, l, u);

    return -1.0;
}
