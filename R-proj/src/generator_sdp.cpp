// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "samplers.h"
#include "rounding.h"
#include "sample_only.h"
#include "sdp_generator.h"
#include "spectrahedron.h"


//' @export
// [[Rcpp::export]]
void generator_sdp(Rcpp::Nullable<int> nn = R_NilValue,
                             Rcpp::Nullable<int> mm = R_NilValue) {

    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Cartesian <NT, NT, VT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef LMI <MT, VT> lmi;
    typedef Spectrahedron <lmi, Point> spectaedro;
    unsigned int n = Rcpp::as<int>(nn);
    //int W = Rcpp::as<unsigned int>(walk_length);

    spectaedro SP;//, SP2;
    SP = generateSDP2<lmi, spectaedro, Point>(Rcpp::as<int>(nn), Rcpp::as<int>(mm));

    //SP2 = SP;

    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    //typedef boost::mt19937 RNGType;
    //RNGType rng(seed);
    //boost::random::uniform_real_distribution<>(urdist);
    //boost::random::uniform_real_distribution<> urdist1(-1, 1);
    //boost::random::uniform_int_distribution<> uidist(0, n - 1);

    //bool round = false;
    //std::pair <Point, NT> InnerB;
    //Point p(Rcpp::as<int>(nn));
    //NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0;
    //InnerB.first = p;// = SP.ComputeInnerBall(diam_spec);

    //vars <NT, RNGType> var(0, Rcpp::as<int>(nn), 1, 1, 0.0, 0.1, 0, 0.0, 0, InnerB.second, diam_spec, rng, urdist,
                           //urdist1,
                           //-1.0, true, false, round, false, false, false, false, false, true);

    //std::list <Point> randPoints, randPoints2;
    //spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
    //settings.LMIatP = SP.getLMI().getA0();
    //preproccess_spectrahedron(SP, p, var, settings, round_value, diam_spec, rad, round);
    //settings.LMIatP = SP.getLMI().getA0();
    //p = Point(n);

    Point c = get_direction<RNGType, Point, NT>(n);

    std::filebuf fb;
    std::string bar = "_";
    std::string txt = ".txt";
    fb.open("sdp_prob"+bar+std::to_string(n)+bar+std::to_string(Rcpp::as<int>(mm))+txt, std::ios::out);
    std::ostream os(&fb);
    writeSDPAFormatFile<MT>(os, SP.getLMI(), c.get_coefficients());

    return;
}

