// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20014-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/ballintersectsimplex.h"
//#include "random_walks/gcw_estimator.hpp"
#include "random_walks/uniform_billiard_gcw_walk.hpp"
#include "random_walks/uniform_gcw_billiard_L.hpp"
#include "diagnostics/psrf_updater.hpp"
#include "diagnostics/univariate_psrf.hpp"



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
Rcpp::NumericMatrix sample_component_psrf_billiard(Rcpp::NumericMatrix A,
                                     Rcpp::NumericVector b,
                                     Rcpp::NumericVector x,
                                     unsigned int N,
                                     unsigned int walk_length,
                                     Rcpp::NumericVector x0,
                                     double psrf_target,
                                     double L)
{
    typedef double NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef UnitBallIntersectSimplex <NT, VT, MT> Body;

    unsigned int d = A.ncol();
    VT b2 = Rcpp::as<VT>(b) - Rcpp::as<MT>(A)*Rcpp::as<VT>(x0);

    

    //MT VV = Rcpp::as<MT>(V);
    //VT Vnorms_shifted(VV.cols());
    //for (int i =0; i<VV.cols(); i++)
    //{
    //    VV.col(i) = VV.col(i) - Rcpp::as<VT>(x0);
    //    Vnorms_shifted(i) = VV.col(i).norm() * VV.col(i).norm();
    //}

    MT V2;// = Rcpp::as<MT>(V);
    //V2 = V2.colwise() - Rcpp::as<VT>(x0);
    //VT Vnorms_shifted = V2.colwise().norm();
    //Vnorms_shifted = Vnorms_shifted.cwiseProduct(Vnorms_shifted);

    //VT Vnorms = Rcpp::as<MT>(V).colwise().norm();
    //Vnorms = Vnorms.cwiseProduct(Vnorms);

    //std::cout<<VV<<"\n\n"<<std::endl;
    //std::cout<<V2<<"\n\n"<<std::endl;

    //std::cout<<Vnorms_shifted<<"\n\n"<<std::endl;
    //std::cout<<Vnorms_shifted2<<"\n\n"<<std::endl;

    Body BS(d, Rcpp::as<MT>(A), b2, V2, VT::Zero(d));
    PSRFestimator<NT, VT, MT> psrf_estimator(N, d);

    //Body BS2(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), Rcpp::as<MT>(V), Rcpp::as<VT>(x0));

    //Body BS3(d, Rcpp::as<MT>(A), Rcpp::as<VT>(b), VV, VT::Zero(d), Vnorms_shifted);

    RNGType rng(d);
    VT p = Rcpp::as<VT>(x), center = Rcpp::as<VT>(x0), y(d);
    p -= center;
    //if (BS2.is_in((p)) == 0)
    //{
    //    std::cout<<"BS2 initial point outside"<<std::endl;
    //    exit(-1);
    //}

    //if (BS.is_in(p) == 0)
   // {
    //    std::cout<<"BS initial point outside"<<std::endl;
    //    exit(-1);
    //}
    //exit(-1);

    typedef BilliardGCWalk::template Walk
            <
                Body,
                RNGType
            > CGWalk;
    
    

    MT samples;//(d, N);
    MT sigma;
    samples.resize(d, N);
    unsigned int iter = 1, MAX_ITER = 1000;
    VT psrf_values(d);
    NT psrf_val;
    
    unsigned int countsIn_total = 0;

    typedef GCWalkL::template Walk
            <
                Body,
                RNGType
            > CGWalkL;
    CGWalkL walkL(BS, p, rng);

    for (int i=0; i<100*d; i++) {
        walkL.template apply(BS, p, walk_length, rng);
    }
    
    NT max_theta = walkL.get_max_theta();
    //std::cout<<"max_theta = "<<max_theta<<std::endl;
    L = 10*max_theta;
    if (L > NT(2)*M_PI) {
        L = NT(2)*M_PI;
    }

    VT p0 = p;

    CGWalk walk(BS, p, L, rng);

    //bool outside = false;
    while (iter <= MAX_ITER)
    {
        countsIn_total = 0;
        p = p0;
        //std::cout<<"initializing..."<<std::endl;
        walk.template initialize(BS, p, rng);
        //std::cout<<"initialized done..."<<std::endl;
        std::cout<<"sampling..."<<std::endl;
        for (int i = 0; i < N; i++)
        {   
            if ((i+1)%1000 == 0) {
                std::cout<<"i = "<<i<<std::endl;
            }
            walk.template apply(BS, p, walk_length, rng);
            samples.col((iter-1)*N + i) = p + center;
            countsIn_total++;

            if (rng.sample_urdist() < (NT(1) / countsIn_total))
            {
                p0 = p;
            }
        }
        //psrf_values = univariate_psrf<NT, VT>(samples);
        //std::cout<<"[1]psrf_values = "<<psrf_values.maxCoeff()<<std::endl;

        psrf_estimator.update_estimator(samples);
        psrf_estimator.estimate_psrf();
        psrf_values = psrf_estimator.get_psrf();
        std::cout<<"[2]psrf_values = "<<psrf_estimator.get_psrf().maxCoeff()<<"\n"<<std::endl;

        psrf_val = psrf_values.maxCoeff();
        if (psrf_val <= psrf_target) {
            //std::cout<<"psrf_val = "<<psrf_val<<std::endl;
            return Rcpp::wrap(samples);
        }
        if (iter == MAX_ITER) {
            return Rcpp::wrap(samples);
        }
        iter++;
        samples.conservativeResize(d, iter*N);
        
    }

    return Rcpp::wrap(samples);
}
