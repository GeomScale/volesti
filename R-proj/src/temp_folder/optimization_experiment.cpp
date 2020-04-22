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
Rcpp::NumericMatrix opti_exp(Rcpp::Nullable<int> nn = R_NilValue,
                               Rcpp::Nullable<int> mm = R_NilValue,
                               Rcpp::Nullable<unsigned int> N = R_NilValue,
                               Rcpp::Nullable<unsigned int> walk_length = R_NilValue){

    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Cartesian<NT, NT, VT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937 RNGType;
    typedef LMI <MT, VT> lmi;
    typedef Spectrahedron<lmi, Point> spectaedro;
    unsigned int n = Rcpp::as<int>(nn);
    int W = Rcpp::as<unsigned int>(walk_length);

    spectaedro SP;//, SP2;
    SP = generateSDP2<lmi, spectaedro, Point>(Rcpp::as<int>(nn), Rcpp::as<int>(mm));

    //SP2 = SP;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    typedef boost::mt19937    RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    boost::random::uniform_int_distribution<> uidist(0, n - 1);

    bool round = false;
    std::pair<Point,NT> InnerB;
    Point p(Rcpp::as<int>(nn));
    NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0;
    InnerB.first = p;// = SP.ComputeInnerBall(diam_spec);

    vars<NT, RNGType> var(0,Rcpp::as<int>(nn), 1, 1,0.0,0.1,0,0.0,0, InnerB.second,diam_spec,rng,urdist,urdist1,
                          -1.0,true,false,round,false,false,false,false,false, true);

    std::list<Point> randPoints, randPoints2;
    spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
    settings.LMIatP = SP.getLMI().getA0();
    preproccess_spectrahedron(SP, p, var, settings, round_value, diam_spec, rad, round);
    settings.LMIatP = SP.getLMI().getA0();
    p = Point(n);

    spectaedro::BoundaryOracleBoltzmannHMCSettings settings2;
    settings2.first = true;
    settings2.epsilon = 0.0001;
    //settings2.LMIatP = SP.getLMI().getA0();
    Point c = get_direction<RNGType, Point, NT>(n);

    std::filebuf fb;
    fb.open ("sdp_prob.txt",std::ios::out);
    std::ostream os(&fb);
    writeSDPAFormatFile<MT>(os, SP.getLMI(), c.get_coefficients());

    NT T = 2.0 * var.diameter, min_val, temp_val;
    std::cout<<"Starting sampling.."<<std::endl;
    VT p_mean(n);
    MT Ratios(4,5);
    for (int k = 0; k < 5; ++k) {

        randPoints.clear();
        randPoints2.clear();
        std::cout << "T = " << T << std::endl;

        //p = Point(n);
        for (int i = 0; i < Rcpp::as<unsigned int>(N); ++i) {
            for (int j = 0; j < Rcpp::as<unsigned int>(walk_length); ++j) {
                HMC_boltzmann_reflections(SP, p, diam_spec, var, c, T, settings2);
            }
            randPoints.push_back(p);
        }
        std::cout << "HMC points sampled.." << std::endl;

        p = Point(n);
        for (int i = 0; i < Rcpp::as<unsigned int>(N); ++i) {
            for (int j = 0; j < Rcpp::as<unsigned int>(walk_length); ++j) {
                hit_and_run_Boltzmann_spec(p, SP, var, c, T);
            }
            randPoints2.push_back(p);
        }
        std::cout << "HNR points sampled.." << std::endl;

        p_mean = VT::Zero(n);
        min_val = std::numeric_limits<NT>::max();
        for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++) {
            p_mean += (*rpit).getCoefficients();
            temp_val = (c.getCoefficients()).dot((*rpit).getCoefficients());
            if (temp_val < min_val) min_val = temp_val;
        }
        p_mean = p_mean * (1.0 / NT(Rcpp::as<unsigned int>(N)));
        Ratios(0,k) = p_mean.dot(c.getCoefficients());
        Ratios(1,k) = min_val;

        p_mean = VT::Zero(n);
        min_val = std::numeric_limits<NT>::max();
        for (typename std::list<Point>::iterator rpit = randPoints2.begin(); rpit!=randPoints2.end(); rpit++) {
            p_mean += (*rpit).getCoefficients();
            temp_val = (c.getCoefficients()).dot((*rpit).getCoefficients());
            if (temp_val < min_val) min_val = temp_val;
        }
        p_mean = p_mean * (1.0 / NT(Rcpp::as<unsigned int>(N)));
        Ratios(2,k) = p_mean.dot(c.getCoefficients());
        Ratios(3,k) = min_val;

        T = T * 0.3;
    }

    return Rcpp::wrap(Ratios);

}