// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
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
Rcpp::NumericMatrix cotrol_exp(Rcpp::Nullable<int> nn = R_NilValue,
                                 Rcpp::Nullable<int> mm = R_NilValue,
                                 Rcpp::Nullable<unsigned int> N = R_NilValue,
                                 Rcpp::Nullable<unsigned int> M = R_NilValue,
                                 Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                                 Rcpp::Nullable<unsigned int> walk_type = R_NilValue){

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

    std::list<Point> randPoints;
    spectaedro::BoundaryOracleBilliardSettings settings(SP.getLMI().getMatricesDim());
    settings.LMIatP = SP.getLMI().getA0();
    preproccess_spectrahedron(SP, p, var, settings, round_value, diam_spec, rad, round);
    settings.LMIatP = SP.getLMI().getA0();
    p = Point(n);

    rand_point_generator_spec(SP, p, 5000, 5, randPoints, var, settings);
    std::cout<<"5000 BW points sampled.."<<std::endl;

    Point c(n);
    MT Cs(3,n), Bs(3,2);

    std::vector<NT> inner_prods;
    std::cout<<"Starting initialization.."<<std::endl;
    for (int i = 0; i < 3; ++i) {
        c = get_direction<RNGType, Point, NT>(n);

        for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++) {
            inner_prods.push_back(c.getCoefficients().dot((*rpit).getCoefficients()));
        }
        std::sort(inner_prods.begin(), inner_prods.end());

        Cs.row(i) = c.getCoefficients();
        Bs(i,0) = inner_prods[int(5000*(0.1+i*0.05))/2.0];
        Bs(i,1) = inner_prods[5000 - int(5000*(0.1+i*0.05))/2.0];

    }
    std::cout<<"Initialization completed.."<<std::endl;
    std::cout<<"Cs = "<<Cs<<"\n\nBs = "<<Bs<<"\n"<<std::endl;

    int walkL = 1, count1, count2, count3, counter = 0;
    VT prods(3);
    int Q = W/5 + 1;
    MT Ratios = MT::Zero(Q*3*3, Rcpp::as<int>(M));
    int jj;

    std::cout<<"Starting BW computations.."<<std::endl;
    while(walkL <= W) {

        std::cout<<"Walk length = "<<walkL<<std::endl;
        for (int i = 0; i <  Rcpp::as<int>(M); ++i) {

            count1 = count2 = count3 = 0;
            randPoints.clear();
            settings.LMIatP = SP.getLMI().getA0();
            p = Point(n);
            rand_point_generator_spec(SP, p, Rcpp::as<int>(N), walkL, randPoints, var, settings);

            jj=0;
            for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) {
                prods = Cs * ((*rpit).getCoefficients());
                if (prods(0) < Bs(0, 0) || prods(0) > Bs(0, 1)) count1++;
                if (jj<90) {
                    if (prods(1) < Bs(1, 0) || prods(1) > Bs(1, 1)) count2++;
                }
                if (jj<110) {
                    if (prods(2) < Bs(2, 0) || prods(2) > Bs(2, 1)) count3++;
                }
            }

            Ratios((counter*3) + 0, i) = NT(count1) / NT(Rcpp::as<int>(N));
            Ratios((counter*3) + 1, i) = NT(count2) / 90.0;
            Ratios((counter*3) + 2, i) = NT(count3) / 110.0;

        }

        if (walkL == 1) {
            walkL = 5;
        } else {
            walkL += 5;
        }
        counter++;
    }

    walkL = 1;
    counter = 0;
    std::cout<<"Starting RDHR computations.."<<std::endl;
    while(walkL <= W) {

        std::cout<<"Walk length = "<<walkL<<std::endl;
        for (int i = 0; i <  Rcpp::as<int>(M); ++i) {

            count1 = count2 = count3 = 0;
            randPoints.clear();
            settings.LMIatP = SP.getLMI().getA0();
            p = Point(n);
            for (int j = 0; j < Rcpp::as<int>(N); ++j) {
                for (int k = 0; k < walkL; ++k) {
                    hit_and_run_spec(p, SP, var);
                }
                randPoints.push_back(p);
            }

            jj=0;
            for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) {
                prods = Cs * ((*rpit).getCoefficients());
                if (prods(0) < Bs(0, 0) || prods(0) > Bs(0, 1)) count1++;
                if (jj<90) {
                    if (prods(1) < Bs(1, 0) || prods(1) > Bs(1, 1)) count2++;
                }
                if (jj<110) {
                    if (prods(2) < Bs(2, 0) || prods(2) > Bs(2, 1)) count3++;
                }
            }

            Ratios((counter*3) + 0 + Q*3, i) = NT(count1) / NT(Rcpp::as<int>(N));
            Ratios((counter*3) + 1 + Q*3, i) = NT(count2) / 90.0;
            Ratios((counter*3) + 2 + Q*3, i) = NT(count3) / 110.0;

        }

        if (walkL == 1) {
            walkL = 5;
        } else {
            walkL += 5;
        }
        counter++;
    }

    walkL = 1;
    counter = 0;
    spectaedro::BoundaryOracleCDHRSettings settings2(SP.getLMI().getMatricesDim());
    settings2.LMIatP = SP.getLMI().getA0();
    NT returnLambda, b;
    Point v(n);
    int rand_coord;
    std::cout<<"Starting CDHR computations.."<<std::endl;
    while(walkL <= W) {

        std::cout<<"Walk length = "<<walkL<<std::endl;
        for (int i = 0; i <  Rcpp::as<int>(M); ++i) {

            count1 = count2 = count3 = 0;
            randPoints.clear();
            settings2.LMIatP = SP.getLMI().getA0();
            settings2.first = true;
            p = Point(n);
            for (int j = 0; j < Rcpp::as<int>(N); ++j) {
                for (int k = 0; k < walkL; ++k) {
                    v = Point(n);
                    rand_coord = uidist(rng);
                    v.set_coord(rand_coord, 1.0);//(rand_coord) = 1.0;
                    std::pair<NT, NT> dbpair = SP.boundaryOracle(p.get_coefficients(), v.get_coefficients());
                    double min_plus = dbpair.first;
                    double max_minus = dbpair.second;
                    Point b1 = (min_plus * v) + p;
                    Point b2 = (max_minus * v) + p;
                    double lambda = urdist(rng);
                    p = (lambda * b1);
                    p = ((1 - lambda) * b2) + p;
                    //hit_and_run_coord_update_spec(p, SP, rand_coord, a, b, var, settings2, returnLambda);
                }
                randPoints.push_back(p);
            }

            jj=0;
            for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++) {
                prods = Cs * ((*rpit).getCoefficients());
                if (prods(0) < Bs(0, 0) || prods(0) > Bs(0, 1)) count1++;
                if (jj<90) {
                    if (prods(1) < Bs(1, 0) || prods(1) > Bs(1, 1)) count2++;
                }
                if (jj<110) {
                    if (prods(2) < Bs(2, 0) || prods(2) > Bs(2, 1)) count3++;
                }
            }

            Ratios((counter*3) + 0 + Q*3*2, i) = NT(count1) / NT(Rcpp::as<int>(N));
            Ratios((counter*3) + 1 + Q*3*2, i) = NT(count2) / 90.0;
            Ratios((counter*3) + 2 + Q*3*2, i) = NT(count3) / 110.0;

        }

        if (walkL == 1) {
            walkL = 5;
        } else {
            walkL += 5;
        }
        counter++;
    }


    return Rcpp::wrap(Ratios);

}