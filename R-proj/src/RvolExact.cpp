// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <random.hpp>
#include <random/uniform_int.hpp>
#include <random/normal_distribution.hpp>
#include <random/uniform_real_distribution.hpp>
#include "polytopes.h"
#include "exact_vols.h"


template <typename FT>
FT factorial(FT n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double Rvol_exact (Rcpp::NumericMatrix A, bool exact_zono, bool exact_cube,
              bool exact_simplex, bool exact_cross, int dim){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;

    typedef Zonotope<Point> zonotope;
    zonotope ZP;

    NT vol;
    if (exact_zono) {
        unsigned int m = A.nrow() - 1;
        unsigned int n = A.ncol() - 1;

        std::vector <std::vector<NT> > Pin(m + 1, std::vector<NT>(n + 1));

        for (unsigned int i=0; i<m+1; i++){
            for(unsigned int j=0; j<n+1; j++){
                Pin[i][j]=A(i,j);
            }
        }

        ZP.init(Pin);
        vol= exact_zonotope_vol<NT>(ZP);
    } else if(exact_simplex) {
        unsigned int m = A.nrow() - 1;
        unsigned int n = A.ncol() - 1;
        if (dim!=n) {
            vol = 1.0 / factorial(NT(dim));
            return vol;
        }
        typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

        MT V(n, n + 1);
        for (int i = 1; i < A.ncol(); ++i) {
            for (int j = 1; j < A.nrow(); ++j) {
                V(i-1, j-1) = A(j, i);
            }
        }
        VT v0 = V.col(n);
        VT V2 = V.block(0,0,n,n);
        V2 = V2.colwise() - v0;

        vol = V2.determinant();
        vol = vol / factorial(NT(n));
    } else if(exact_cube) {
        vol = std::pow(2.0, NT(dim));
    } else {
        vol = std::pow(2.0, NT(dim));
        vol = vol / factorial(NT(dim));
    }

    return vol;

}