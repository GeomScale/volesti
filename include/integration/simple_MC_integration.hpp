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

enum volumeType { CB , CG , SOB }; // Volume type for polytope
typedef const unsigned int Uint; // positive constant value for no of samples & dimensions
VT Origin(0); // undefined 0 dimensional VT

// To check if two n-dimensional points ensure valid limits in integration 
bool validLimits(VT LL, VT UL){
    if( UL.rows() == LL.rows() ) {
        for ( int i = 0 ; i< LL.rows() ; i++){
            if( LL(i) > UL(i) ){
            std::cout << "Invalid integration limits\n";
            return false;
            }
        }
        return true;
    }else{
        std::cout << "Invalid integration limits\n";
        return false;
    }
}

// Hyper-Rectangle volume calculator in n-dimensions
NT hyperRectVolume(VT LL, VT UL){
    NT product=1;
    if(validLimits(LL,UL)){
        for(int i=0; i<LL.rows(); ++i){
            product = product * ( UL(i) - LL(i) );
        }
        return product;
    }
    else return 0;
}

// To sample a point between two n-dimensional points
VT samplerBWLimits(VT LL, VT UL){
    VT sample_point(LL.rows()); NT x;
    for(int i=0; i<LL.rows(); ++i){
        sample_point(i) = LL(i) + (NT)(rand()) / ((NT)(RAND_MAX/(UL(i) - LL(i))));     
    }
    return sample_point;
}

// Simple MC Integration over Hyper-Rectangles
template < typename Functor >
void simple_mc_integrate(Functor Fx, Uint N ,VT LL, VT UL){
    NT sum = 0;
    if(validLimits(LL,UL)){
        for (int i = 1; i < N; i++ ) {
            sum = sum + Fx(samplerBWLimits(LL,UL));
        }    
        std::cout << "Integral Value: " << hyperRectVolume(LL,UL) * sum / N << "\n";
    }else{
        std::cout << "Invalid integration limits\n";
    }
}

template < typename Functor , typename Polytope >
void simple_mc_polytope_integrate(Functor Fx, Polytope &P, Uint N,volumeType vType=CB,VT newOrigin=Origin){

    // Polytope volumetric calculation params
    Uint dim = P.dimension();
    int walk_length = 10 + dim/10; 
    NT e=0.1;
    NT volume=0;

    // Volume calculation for HPolytope
    switch(vType){
        case CB:     
            volume = volume_cooling_balls<BallWalk, RNGType, HPOLYTOPE>(P, e, walk_length).second;
            break;
        case CG: 
            volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(P, e, walk_length);
            break;
        case SOB: 
            volume = volume_sequence_of_balls<BallWalk, RNGType>(P, e, walk_length);
            break;    
        default:  
            volume = 0;
            std::cout << "Error in enum:vType\n";
            break;
    }
    
    // For implementation of ReHMC walks
    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();
    RNGType rng(1);
    Point x0 = inner_ball.first;
    NT r = inner_ball.second;
    GaussianFunctor::parameters<NT, Point> params(x0, 2 / (r * r), NT(-1));
    GaussianRDHRWalk::Walk<HPOLYTOPE, RNGType> walk(P, x0, params.L, rng);
    NegativeGradientFunctor F(params);
    NegativeLogprobFunctor f(params);
    HamiltonianMonteCarloWalk::parameters<NT, NegativeGradientFunctor> hmc_params(F, dim);
    HamiltonianMonteCarloWalk::Walk
    <Point,HPOLYTOPE, RNGType, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
    hmc( &P , x0, F, f, hmc_params);

    // Check if origin is shifted
    const bool &shiftedOrigin = ( newOrigin.rows()==0 )? false : true ;

    // Taking samples using ReHMC walks
    MT samples(N,dim);
    if(shiftedOrigin){
        if(newOrigin.rows() == dim){
            for (int i = 0; i < N; i++ ) {
                hmc.apply(rng, walk_length);
                samples.row(i) = hmc.x.getCoefficients() + newOrigin;
            }
        }else{
            volume = 0;
            std::cout << "Error in Polytope & shiftedOrigin dimensions do not match\n";
        }
    }else{
        for (int i = 0; i < N; i++ ) {
            hmc.apply(rng, walk_length);
            samples.row(i) = hmc.x.getCoefficients();
        }
    }

    // Evaluation of sampled points
    NT sum = 0;
    for (int i = 0; i < N; i++ ) {
        sum = sum + Fx(samples.row(i));
    }

    // Final step for integration
    std::cout << "Integral Value over H-Polytope: " << volume* sum / N << "\n";    

}

#endif
