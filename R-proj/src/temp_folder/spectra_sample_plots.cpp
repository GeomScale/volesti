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
Rcpp::NumericMatrix spectra_plot(Rcpp::Nullable<int> nn = R_NilValue,
                                  Rcpp::Nullable<int> mm = R_NilValue,
                                  Rcpp::Nullable<unsigned int> N1 = R_NilValue,
                                  Rcpp::Nullable<unsigned int> N2 = R_NilValue,
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

    spectaedro SP;//, SP2;
    SP = generateSDP2<lmi, spectaedro, Point>(Rcpp::as<int>(nn), Rcpp::as<int>(mm));
    //SP2 = SP;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    typedef boost::mt19937    RNGType;
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    bool round = false;
    std::pair<Point,NT> InnerB;
    Point p(Rcpp::as<int>(nn));
    NT nballs2, diam_spec, vol_spec, rad, round_value = 1.0;
    InnerB.first = p;// = SP.ComputeInnerBall(diam_spec);

    vars<NT, RNGType> var(0,Rcpp::as<int>(nn), 1, 1,0.0,0.1,0,0.0,0, InnerB.second,diam_spec,rng,urdist,urdist1,
                          -1.0,true,false,round,false,false,false,false,false, true);

    std::list<Point> randPoints, randPoints2, randPoints3;
    spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
    settings.LMIatP = SP.getLMI().getA0();

    preproccess_spectrahedron(SP, p, var, settings, round_value, diam_spec, rad, round);
    p = Point(n);
    boundary_rand_point_generator_spec(SP, p, Rcpp::as<unsigned int>(N1), Rcpp::as<unsigned int>(walk_length), randPoints, var);
    std::cout<<"Boundary points sampled.."<<std::endl;
    
    settings.LMIatP = SP.getLMI().getA0();
    p = Point(n);
    rand_point_generator_spec(SP, p, Rcpp::as<unsigned int>(N2), Rcpp::as<unsigned int>(walk_length), randPoints2, var, settings);
    std::cout<<"Interior points sampled.."<<std::endl;


    spectaedro::BoundaryOracleBoltzmannHMCSettings settings2;
    settings2.first = true;
    settings2.epsilon = 0.0001;
    //settings2.LMIatP = SP.getLMI().getA0();
    p = Point(n);
    Point c = get_direction<RNGType, Point, NT>(n);
    std::cout<<c.get_coefficients()<<std::endl;
    NT T = 1.0;
    std::cout<<"Starting HMC sampling.."<<std::endl;
    for (int i = 0; i < Rcpp::as<unsigned int>(N2); ++i) {
        for (int j = 0; j < Rcpp::as<unsigned int>(walk_length); ++j) {
            HMC_boltzmann_reflections(SP, p, diam_spec, var, c, T, settings2);
        }
        randPoints3.push_back(p);
    }
    std::cout<<"HMC points sampled..\n"<<std::endl;

    std::filebuf fb;
    fb.open ("sdp_prob.txt",std::ios::out);
    std::ostream os(&fb);
    writeSDPAFormatFile<MT>(os, SP.getLMI(), c.get_coefficients());
    //SP.print_matrices();

    MT RetMat(n, 2*Rcpp::as<unsigned int>(N1)+ 2*Rcpp::as<unsigned int>(N2));
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
        RetMat.col(jj) = (*rpit).getCoefficients();

    jj = 0;
    for (typename std::list<Point>::iterator rpit = randPoints2.begin(); rpit!=randPoints2.end(); rpit++, jj++)
        RetMat.col(jj+2*Rcpp::as<unsigned int>(N1)) = (*rpit).getCoefficients();

    jj = 0;
    for (typename std::list<Point>::iterator rpit = randPoints3.begin(); rpit!=randPoints3.end(); rpit++, jj++)
        RetMat.col(jj+2*Rcpp::as<unsigned int>(N1)+Rcpp::as<unsigned int>(N2)) = (*rpit).getCoefficients();

    return Rcpp::wrap(RetMat);

}