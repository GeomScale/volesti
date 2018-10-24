// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
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
#include "gaussian_samplers.h"
#include "gaussian_annealing.h"
#include "zonovol_heads/sampleTruncated.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix vol_zono (Rcpp::Reference P, double e, Rcpp::Function mvrandn) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef VPolytope <Point, RNGType> Vpolytope;
    typedef Zonotope <Point> zonotope;
    //unsigned int n_threads=1,i,j;

    bool rand_only = false,
            NN = false,
            birk = false,
            verbose = false;
    unsigned int n_threads = 1;
    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;

    int N = 1000;
    Rcpp::NumericVector l = Rcpp::NumericVector::create(-1.0, -1.0);
    Rcpp::NumericVector u = Rcpp::NumericVector::create(1.0, 1.0);
    Rcpp::NumericMatrix sig(2, 2);
    sig(0, 0) = 2, sig(0, 1) = 0;
    sig(1, 0) = 0;
    sig(1, 1) = 2;
    //Rcpp::NumericMatrix::Column zzcol = sig( _, 1);
    //Point p(2,std::vector<NT>::iterator(sig( _, 0).begin()), std::vector<NT>::iterator(sig( _, 0).end()));
    //p.print();
    //std::cout<<p[0]<<" "<<"dfsdf "<<p[1]<<std::endl;

            std::list<Point> randPoints;
    randPoints = sampleTr(l, u, sig, mvrandn, randPoints);
    typename std::list<Point>::iterator pit=randPoints.begin();
    for ( ;  pit!=randPoints.end(); ++pit) {
        (*pit).print();
    }
    //std::vector<std::vector<NT> > Pin(X.begin(), X.end());
    return sig;
}/*
    Rcpp::NumericMatrix A = P.field("G");

    unsigned int m=A.nrow()-1;
    unsigned int n=A.ncol()-1;

    std::vector<std::vector<NT> > Pin(m+1, std::vector<NT>(n+1));

    for (unsigned int i=0; i<m+1; i++){
        for(unsigned int j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
        }
    }

    ZP.init(Pin);
    //Eigen::MatrixXd pinv = A.completeOrthogonalDecomposition().pseudoInverse();
    std::pair<Point,NT> InnerB;
    InnerB = ZP.ComputeInnerBall();

}*/