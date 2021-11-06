// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20014-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/ballintersectsimplex.h"
#include "random_walks/gcw_estimator.hpp"
#include "random_walks/uniform_great_cycle_walk.hpp"

//' Gelman-Rubin and Brooks-Gelman Potential Scale Reduction Factor (PSRF) for each marginal
//'
//' @param samples A matrix that contans column-wise the sampled points from a geometric random walk.
//' @param method A string to reauest diagnostic: (i) \code{'normal'} for psrf of Gelman-Rubin and (ii) \code{'interval'} for psrf of Brooks-Gelman.
//'
//' @references \cite{Gelman, A. and Rubin, D. B.,
//' \dQuote{Inference from iterative simulation using multiple sequences,} \emph{Statistical Science,} 1992.}
//'
//' @references \cite{Brooks, S. and Gelman, A.,
//' \dQuote{General Methods for Monitoring Convergence of Iterative Simulations,} \emph{Journal of Computational and Graphical Statistics,} 1998.}
//'
//' @return A vector that contains the values of PSRF for each coordinate
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_component(Rcpp::NumericMatrix A,
                                     Rcpp::NumericVector b,
                                     Rcpp::NumericVector x,
                                     unsigned int N,
                                     unsigned int walk_length,
                                     Rcpp::NumericVector x0,
                                     Rcpp::NumericMatrix V,
                                     Rcpp::NumericVector Vnorms,
                                     double c1,
                                     double c2)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef UnitBallIntersectSimplex <NT, VT, MT> Body;

    unsigned int d = A.ncol();
    VT b2 = Rcpp::as<VT>(b) + Rcpp::as<MT>(A)*Rcpp::as<VT>(x0);

    Body BS(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), Rcpp::as<MT>(V), Rcpp::as<VT>(x0), Rcpp::as<VT>(Vnorms));

    Body BS2(d, Rcpp::as<MT>(A), b2, Rcpp::as<MT>(V), Rcpp::as<VT>(x0), Rcpp::as<VT>(Vnorms));

    MT VV = Rcpp::as<MT>(V);
    VT Vnorms_shifted(VV.cols());
    for (int i =0; i<VV.cols(); i++)
    {
        VV.col(i) = VV.col(i) - Rcpp::as<VT>(x0);
        Vnorms_shifted(i) = VV.col(i).norm() * VV.col(i).norm();
    }
    std::cout<<VV<<"\n"<<std::endl;

    Body BS3(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), VV, VT::Zero(d), Vnorms_shifted);

    RNGType rng(d);
    VT p = Rcpp::as<VT>(x), center = Rcpp::as<VT>(x0), y(d);

    if (BS2.is_in((p+center)) == 0)
    {
        std::cout<<"BS2 initial point outside"<<std::endl;
        exit(-1);
    }

    if (BS3.is_in(p) == 0)
    {
        std::cout<<"BS3 initial point outside"<<std::endl;
        exit(-1);
    }
    //exit(-1);

    typedef GCWalk::template Walk
            <
                Body,
                RNGType
            > CGWalk;
    
    CGWalk walk(BS, p, rng);

    MT samples(d, N);


    for (int i = 0; i < N; i++)
    {
        walk.template apply(BS, p, walk_length, rng);
        y = p + center;
        if (BS2.is_in(y) == 0)
        {
            std::cout<<"BS2 point outside"<<std::endl;
            exit(-1);
        }
        if (BS3.is_in(p) == 0)
        {
            std::cout<<"BS3 point outside"<<std::endl;
            exit(-1);
        }
        samples.col(i) = p + center;
    }
    typedef std::vector<VT> PointList;

    typedef GCWEstimator::template Walk
            <
                Body,
                PointList,
                RNGType
            > CGEstimator;
    
    CGEstimator estimator(BS3, p, rng);
    estimator.activate_storing();
    NT val = 0, Ntot = 1000, ratio = 1;
    unsigned int W = 1000;
    
    PointList list_of_points = estimator.estimate(BS3,
                       BS3,
                       p,
                       walk_length,
                       0.1,
                       val,
                       W,
                       Ntot,
                       ratio,
                       rng);
    std::cout<<"val = "<<val<<std::endl;
    std::cout<<"list_of_points.size() = "<<list_of_points.size()<<std::endl;
    int counter1 = 0, counter2 = 0;
    for (int i =0; i<list_of_points.size(); i++)
    {
        y = list_of_points[i] + center;
        if (BS2.is_in(y) == 0)
        {
            std::cout<<"BS2 point outside"<<std::endl;
            exit(-1);
        }
        else{
            counter1++;
        }
        if (BS3.is_in(list_of_points[i]) == 0)
        {
            std::cout<<"BS3 point outside"<<std::endl;
            exit(-1);
        }
        else{
            counter2++;
        }
    }
    std::cout<<"counter1 = "<<counter1<<", counter2 = "<<counter2<<std::endl;

    return Rcpp::wrap(samples);
}
