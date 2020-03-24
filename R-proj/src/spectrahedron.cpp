//
// Created by panagiotis on 2/23/20.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "LMI.h"
#include "spectrahedron.h"
#include "SDPAFormatManager.h"



//' @export
// [[Rcpp::export]]
void writeSdpaFile(Rcpp::Nullable<Rcpp::Reference> spectrahedron = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> objectiveFunction = R_NilValue,
               Rcpp::Nullable<std::string> outputFile = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI <NT, MT, VT> LMI;
    typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;

    std::vector<MT> matrices = Rcpp::as<std::vector<MT> > (Rcpp::as<Rcpp::Reference> (spectrahedron).field("matrices"));
    LMI lmi(matrices);
    SPECTRAHEDRON _spectrahedron(lmi);
    Point c(Rcpp::as<VT> (objectiveFunction));

    std::filebuf fb;
    fb.open(Rcpp::as<std::string> (outputFile), std::ios::out);
    std::ostream os(&fb);

    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.writeSDPAFormatFile(os, _spectrahedron, c);

    return;
}