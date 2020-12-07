// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#include <Rcpp.h>
#include <RcppEigen.h>
//#include <chrono>
//#include "cartesian_geom/cartesian_kernel.h"
//#include <boost/random.hpp>
//#include <boost/random/uniform_int.hpp>
//#include <boost/random/normal_distribution.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include "random_walks/random_walks.hpp"
//#include "volume/volume_sequence_of_balls.hpp"
//#include "volume/volume_cooling_gaussians.hpp"
//#include "preprocess/get_full_dimensional_polytope.hpp"
//#include "extractMatPoly.h"


// [[Rcpp::export]]
Rcpp::NumericVector solve_undetermined_system_lu (Rcpp::NumericMatrix Ar, Rcpp::NumericVector br)
{
    typedef double NT;
    //typedef Cartesian<NT>    Kernel;
    //typedef typename Kernel::Point    Point;
    //typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    MT Aeq = Rcpp::as<MT>(Ar);
    VT beq = Rcpp::as<VT>(br);

    VT p = Aeq.fullPivLu().solve(beq);
    
    return Rcpp::wrap(p);
}
