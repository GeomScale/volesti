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

//' Write a SDPA format file
//'
//' Outputs a spectrahedron (the matrices defining a linear matrix inequality) and a vector (the objective function)
//' to a SDPA format file.
//'
//' @param spectrahedron A spectrahedron in n dimensions; must be an object of class Spectrahedron
//' @param objective_function A numerical vector of length n
//' @param output_file Name of the output file
//'
//' @examples
//' \donttest{
//' A0 = matrix(c(-1,0,0,0,-2,1,0,1,-2), nrow=3, ncol=3, byrow = TRUE)
//' A1 = matrix(c(-1,0,0,0,0,1,0,1,0), nrow=3, ncol=3, byrow = TRUE)
//' A2 = matrix(c(0,0,-1,0,0,0,-1,0,0), nrow=3, ncol=3, byrow = TRUE)
//' lmi = list(A0, A1, A2)
//' S = Spectrahedron(matrices = lmi)
//' objFunction = c(1,1)
//' write_sdpa_format_file(S, objFunction, "output.txt")
//' }
//' @export
// [[Rcpp::export]]
void write_sdpa_format_file(Rcpp::Reference spectrahedron,
                            Rcpp::NumericVector objective_function,
                            std::string output_file) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI<NT, MT, VT> LMI;
    typedef Spectrahedron<NT, MT, VT> Spectrahedron;

    std::vector<MT> matrices = Rcpp::as<std::vector<MT> > (spectrahedron.slot("matrices"));
    LMI lmi(matrices);
    Spectrahedron _spectrahedron(lmi);
    Point c(Rcpp::as<VT> (objective_function));

    std::filebuf fb;
    fb.open(output_file, std::ios::out);
    std::ostream os(&fb);

    SdpaFormatManager<NT> sdpaFormatManager;
    sdpaFormatManager.writeSDPAFormatFile(os, _spectrahedron, c);

    return;
}

