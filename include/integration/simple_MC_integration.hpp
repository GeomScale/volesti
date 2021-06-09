// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SIMPLE_MC_INTEGRATION_HPP
#define SIMPLE_MC_INTEGRATION_HPP

#include "hpolytope.h"
#include "vpolytope.h"
#include "Eigen/Eigen"
#include "known_polytope_generators.h"
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <boost_random_number_generator.hpp>
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "misc.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef std::vector<Point> pts;
typedef HPolytope<Point> HPOLYTOPE;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef GaussianFunctor::GradientFunctor<Point> NegativeGradientFunctor;
typedef GaussianFunctor::FunctionFunctor<Point> NegativeLogprobFunctor;
typedef LeapfrogODESolver<Point, NT, HPOLYTOPE , NegativeGradientFunctor> Solver;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT;

typedef const unsigned int Uint; // positive constant value for no of samples & dimensions
typedef std::vector<NT> Vect;

// To return ||X||^2 for a VT
NT normSquared(VT vt){
    NT sum=0;
    for( int i=0; i < vt.cols() ; i++ ){
        sum = sum + vt(i)*vt(i);
    }
    return sum;
}

// To check if two n-dimensional points ensure proper limits in integration 
bool legitLimits(VT LL, VT UL){
   if( UL.cols() == LL.cols() ) {
       for ( int i = 0 ; i< LL.size() ; i++){
           if( LL(i) > UL(i) ){
               std::cout << "Enter limits correctly!\n";
               return false;
            }
       }
       return true;
   }else{
       std::cout << "Enter limits correctly!\n";
       return false;
   }
}

// Hyper-Rectangle volume calculator in n-dimensions
NT hyperRectVolume(VT LL, VT UL){
    NT product=1;
    if(legitLimits(LL,UL)){
        for(int i=0; i<LL.size(); ++i){
            product = product * abs(UL(i) - LL(i));
        }
        return product;
    }
    else return 0;
}

// To sample a point between two n-dimensional points
VT samplerBWLimits(VT LL, VT UL){
    VT sample_point(LL.cols()); NT x;
    for(int i=0; i<LL.cols(); ++i){
        sample_point(i) = LL(i) + (NT)(rand()) / ((NT)(RAND_MAX/(UL(i) - LL(i))));      
    }
    return sample_point;
}

// Simple MC Integration over Hyper-Rectangles
template < typename Functor >
void SimpleMCIntegrate(Functor Fx, Uint N ,VT LL, VT UL){
    NT sum = 0;
    if(legitLimits(LL,UL)){
        for (int i = 1; i < N; i++ ) {
            sum = sum + Fx(samplerBWLimits(LL,UL));
        }    
        std::cout << "Integral Value: " << hyperRectVolume(LL,UL) * sum / N << "\n";
    }else{
        std::cout << "Enter integration limits properly\n";
    }
}

VT Origin(0);
template < typename Functor >
void SimpleMCPolytopeIntegrate(Functor Fx, Uint dim, Uint N ,VT newOrigin=Origin){

    // Polytope volumetric calculation params
    int walk_length = 10 + dim/10; 
    NT e=0.1;
    NT volume=0;

    // Creating a HPolytope and calculating the volume
    HPOLYTOPE HP = generate_cube<HPOLYTOPE>(dim, false);
    //HPOLYTOPE HP = generate_simplex<HPOLYTOPE>(dimensions, false);

    // Volume calculation for HPolytope
    volume = volume_cooling_balls<BallWalk, RNGType, HPOLYTOPE>(HP, e, walk_length).second;
    // volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, e, walk_len);
    // volume = volume_sequence_of_balls<BallWalk, RNGType>(HP, e, walk_len);
    
    // For implementation of ReHMC walks
    std::pair<Point, NT> inner_ball = HP.ComputeInnerBall();
    RNGType rng(1);
    Point x0 = inner_ball.first;
    NT r = inner_ball.second;
    GaussianFunctor::parameters<NT, Point> params(x0, 2 / (r * r), NT(-1));
    GaussianRDHRWalk::Walk<HPOLYTOPE, RNGType> walk(HP, x0, params.L, rng);
    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);
    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);
    HamiltonianMonteCarloWalk::Walk
    <Point,HPOLYTOPE, RNGType, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
    hmc( &HP , x0, F, f, hmc_params);

    // Taking samples using ReHMC walks
    MT samples(N,dim);
    for (int i = 0; i < N; i++ ) {
        hmc.apply(rng, walk_length);
        samples.row(i) = hmc.x.getCoefficients();
    }

    // Check if origin is shifted
    const bool &shiftedOrigin = (newOrigin.rows()==0 )? false : true ;

    // Evaluation of sampled points
    NT sum = 0;
    if(shiftedOrigin){
        if(newOrigin.rows() == dim){
            for (int i = 0; i < N; i++ ) {
            sum = sum + Fx(samples.row(i)+newOrigin);
            }
        }else{
            std::cout << "Dimensions should match to the new origin!"
        }
    }else{
        for (int i = 0; i < N; i++ ) {
        sum = sum + Fx(samples.row(i));
        }
    }

    // Final step for integration
    std::cout << "Integral Value over H-Polytope: " << volume * sum / N << "\n";    

}

#endif