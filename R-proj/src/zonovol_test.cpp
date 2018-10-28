

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
#define VOLESTI_DEBUG
#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(BH)]]
#define VOLESTI_DEBUG
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
#include "zonovol_heads/outer_zono.h"
#include "zonovol_heads/cg_zonovol.h"



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double vol_zono (Rcpp::Reference P, double e, Rcpp::Function mvrandn, bool verbose, double delta_in=0.0, double var_in=0.0, double up_lim=0.3) {

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
    zonotope ZP;
    Rcpp::NumericMatrix A = P.field("G");

    unsigned int m=A.nrow();
    unsigned int n=A.ncol();
    //std::cout<<"hello"<<std::endl;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    MT V = Rcpp::as<MT>(A);
    VT vec = VT::Ones(m);
    ZP.init(n, V, vec);
    MT G = V.transpose();
    MT ps = G.completeOrthogonalDecomposition().pseudoInverse();
    MT sigma = ps*ps.transpose();
    //std::cout<<sigma<<std::endl;
    for (int i1 = 0; i1 < m; ++i1) {
        sigma(i1,i1) = sigma(i1,i1) + 0.00000001;
    }
    //sigma = sigma + 0.00001*MT::DiagonalMatrix(m);

    VT l = VT::Zero(m) - VT::Ones(m);
    VT u = VT::Ones(m);

    //MT test = sampleTr(l, u, sigma, 10, mvrandn, G);
    std::list<Point> randPoints;
    NT delta, ratio;
    double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
    get_delta<Point>(ZP, l, u, sigma, mvrandn, G, var_in, delta_in, up_lim, ratio, randPoints);
    double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(verbose) std::cout << "[1] computation of outer time = " << tstop1 - tstart1 << std::endl;
    delta = delta_in;

    std::pair<Point,NT> InnerB = ZP.ComputeInnerBall();
    NT C = 2.0, frac=0.1, ratio2 = 1.0-1.0/(NT(n));
    int W = 4*n*n+500, N = 500 * ((int) C) + ((int) (n * n / 2));

    vars<NT, RNGType> var2(1, n, 10 + n / 10, 1, 0.0, e, 0, 0.0, 0, InnerB.second, rng,
                           urdist, urdist1, -1.0, verbose, false, false, NN, birk, false,
                           false);
    vars_g<NT, RNGType> var1(n, 1, N, W, 1, e/2.0, InnerB.second, rng, C, frac,
                             ratio2, -1.0, false, verbose, false, false, NN, birk,
                             false, false);

    //sigma = variance * sigma;
    l = l - VT::Ones(m) * delta;
    u = u + VT::Ones(m) * delta;
    if (verbose) std::cout<<"delta = "<<delta<<" first estimation of ratio (t-test value) = "<<ratio<<std::endl;
    //std::cout<<l<<"\n"<<u<<std::endl;
    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
    NT vol = cg_volume_zono(ZP, var1, var2, InnerB, mvrandn, sigma, l, u);
    double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(verbose) std::cout << "[2] outer volume estimation with cg algo time = " << tstop2 - tstart2 << std::endl;
    if (verbose) std::cout<<"volume of outer = "<<vol<<std::endl;

    NT countIn = ratio*1200.0, totCount = 1200.0;
    //std::cout<<"countIn = "<<countIn<<std::endl;

    sigma = 2*var_in * sigma;
    //std::cout<<sigma<<std::endl;
    int count;
    double tstart3 = (double)clock()/(double)CLOCKS_PER_SEC;
    //MT sample = sampleTr(l, u, sigma, 8800, mvrandn, G, count);
    //countIn = countIn + NT(8800-count);
    //totCount = totCount + NT(8800-count);
    MT sample = sampleTr(l, u, sigma, 8800, mvrandn, G);
    randPoints.clear();
    //for (int i = 0; i < count; ++i) {
    for (int i = 0; i < 8800; ++i) {
        Point p(n, typename std::vector<NT>::iterator(sample.col(i).data()), typename std::vector<NT>::iterator(sample.col(i).data() + n));
        randPoints.push_back(p);
    }
    std::list<Point>::iterator rpit = randPoints.begin();

    for ( ;  rpit!=randPoints.end(); ++rpit) {
        if(ZP.is_in(*rpit)==-1) {
            countIn = countIn + 1.0;
        }
        totCount = totCount + 1.0;
    }
    if (verbose) std::cout<<"countIn = "<<countIn<<" totCountIn = "<<totCount<<std::endl;
    if (verbose) std::cout<<"variance = "<<var_in<<std::endl;
    double tstop3 = (double)clock()/(double)CLOCKS_PER_SEC;
    if(verbose) std::cout << "[3] rejection time = " << tstop3 - tstart3 << std::endl;
    if (verbose) std::cout<<"final ratio = "<<countIn / totCount<<std::endl;
    vol = vol * (countIn / totCount);
    //std::cout<<"final volume = "<<vol<<std::endl;
    return vol;
}
