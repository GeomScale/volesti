// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "convex_bodies/spectrahedra/LMI.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"

//'  An internal Rccp function to read a SDPA format file
//'
//' @param input_file Name of the input file
//'
//' @keywords internal
//'
//' @return A list with two named items: an item "matrices" which is a list of the matrices and an vector "objFunction"
// [[Rcpp::export]]
Rcpp::List load_sdpa_format_file(Rcpp::Nullable<std::string> input_file = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI<NT, MT, VT> LMI;
    typedef Spectrahedron<NT, MT, VT> Spectrahedron;

    Spectrahedron _spectrahedron;
    Point c;

    // open stream
    std::ifstream os;
    os.open(Rcpp::as<std::string> (input_file),std::ifstream::in);

    // read file
    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.loadSDPAFormatFile(os, _spectrahedron, c);

    std::vector<MT> matrices = _spectrahedron.getLMI().getMatrices();

    // return spectrahedron and objective function
    Rcpp::List _matrices;

    for (auto matrix : matrices)
        _matrices.push_back(Rcpp::wrap(matrix));

    Rcpp::List retList = Rcpp::List::create(Rcpp::Named("matrices") = _matrices , Rcpp::_["objFunction"] = Rcpp::wrap(c.getCoefficients()));
    return Rcpp::wrap(retList);
}

