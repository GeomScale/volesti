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
Rcpp::List estimate_component(Rcpp::NumericMatrix A,
                                       Rcpp::NumericVector b,
                                       Rcpp::NumericVector x,
                                       unsigned int walk_length,
                                       unsigned int win_len,
                                       Rcpp::NumericVector x0,
                                       Rcpp::NumericMatrix V,
                                       Rcpp::NumericVector c11,
                                       Rcpp::NumericVector c22,
                                       double error,
                                       double ratio,
                                       unsigned int Ntot,
                                       bool storing)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef UnitBallIntersectSimplex <NT, VT, MT> Body;

    unsigned int d = A.ncol();
    RNGType rng(d);

    VT c1 = Rcpp::as<VT>(c11), c2 = Rcpp::as<VT>(c22);

    VT b1 = c1.cwiseProduct(Rcpp::as<VT>(b)) - Rcpp::as<MT>(A)*Rcpp::as<VT>(x0);
    VT b2 = c2.cwiseProduct(Rcpp::as<VT>(b)) - Rcpp::as<MT>(A)*Rcpp::as<VT>(x0);

    //MT V2 = c2 * Rcpp::as<MT>(V);
    //V2 = V2.colwise() - Rcpp::as<VT>(x0);
    //VT Vnorms_shifted = V2.colwise().norm();
    //Vnorms_shifted = Vnorms_shifted.cwiseProduct(Vnorms_shifted);

    //for (int i =0; i<V2.cols(); i++)
    //{
    //    V2.col(i) = V2.col(i) - Rcpp::as<VT>(x0);
    //    Vnorms_shifted(i) = V.col(i).dot(V.col(i));
    //}

    //std::cout<<V2<<"\n\n"<<std::endl;
    //std::cout<<(c2 * Rcpp::as<MT>(V)).colwise() - Rcpp::as<VT>(x0)<<"\n\n"<<std::endl;

    Body BS1(d, Rcpp::as<MT>(A), b1, (c1 * Rcpp::as<MT>(V)).colwise() - Rcpp::as<VT>(x0), VT::Zero(d));
    Body BS2(d, Rcpp::as<MT>(A), b2, (c2 * Rcpp::as<MT>(V)).colwise() - Rcpp::as<VT>(x0), VT::Zero(d));

    VT p = Rcpp::as<VT>(x), center = Rcpp::as<VT>(x0), y(d);
    p -= center;
   
    typedef std::vector<VT> PointList;
    typedef GCWEstimator::template Walk
            <
                Body,
                PointList,
                RNGType
            > CGEstimator;
    //std::cout<<"[1] BS2 point outside"<<std::endl;
    CGEstimator estimator(BS1, p, rng);
    if(storing)
    {
        estimator.activate_storing();
    }
    NT val = 0;
    //std::cout<<"BS2 point outside"<<std::endl;
    PointList list_of_points = estimator.estimate(BS1,
                       BS2,
                       p,
                       walk_length,
                       error,
                       val,
                       win_len,
                       Ntot,
                       ratio,
                       rng);
    //std::cout<<"val = "<<val<<std::endl;
    //std::cout<<"list_of_points.size() = "<<list_of_points.size()<<std::endl;
    //int counter1 = 0, counter2 = 0;
    MT samples;
    if (storing)
    {
        samples.setZero(d, list_of_points.size());
        for (int i =0; i<list_of_points.size(); i++)
        {
            //y = list_of_points[i];
            //if (BS2.is_in(y) == 0)
            //{
            //    std::cout<<"BS2 point outside"<<std::endl;
            //    exit(-1);
            //}
            //else{
                //counter1++;
            samples.col(i) = list_of_points[i] + center;
            //}
        }
    }
    //std::cout<<"counter1 = "<<counter1<<", counter2 = "<<counter2<<std::endl;

    return Rcpp::List::create(Rcpp::Named("ratio") = val, Rcpp::Named("samples") = Rcpp::wrap(samples));
}
