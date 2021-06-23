// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed and/or modified by Suraj Choubey, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SIMPLE_MC_INTEGRATION_HPP
#define SIMPLE_MC_INTEGRATION_HPP

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "convex_bodies/hpolytope.h"
#include "Eigen/Eigen"
#include "generators/known_polytope_generators.h"
#include "boost_random_number_generator.hpp"
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
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT; 

typedef const unsigned int Uint;  // positive constant value for no of samples & dimensions
enum volumeType { CB , CG , SOB }; // Volume type for polytope

// To check if two n-dimensional points ensure valid limits in integration 
bool validLimits(Point LL, Point UL){
    if( UL.dimension() == LL.dimension() ) {
        for ( int i = 0 ; i< LL.dimension() ; i++){
            if( LL[i] > UL[i]){
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
NT hyperRectVolume(Point LL, Point UL){
    NT product=1;
    if(validLimits(LL,UL)){
        for(int i=0; i<LL.dimension(); ++i){
            product = product * ( UL[i] - LL[i] );
        }
        return product;
    }
    else return 0;
}

// To sample a point between two n-dimensional points
Point samplerBWLimits(Point LL, Point UL){
    Point sample_point(LL.dimension()); NT x;
    for(int i=0; i<LL.dimension(); ++i){
        sample_point.set_coord(i , LL[i] + (NT)(rand()) / ((NT)(RAND_MAX/(UL[i] - LL[i]))) );
    }
    return sample_point;
}

// Simple MC Integration over Hyper-Rectangles
template
<
    typename Functor
>
void simple_mc_integrate(Functor Fx, Uint N ,Point LL, Point UL){
    NT sum = 0;
    if(validLimits(LL,UL)){
        for (int i = 0; i <=N; i++ ) {
            sum = sum + Fx(samplerBWLimits(LL,UL));
        }    
        std::cout << "Integral Value: " << hyperRectVolume(LL,UL) * sum / N << "\n"; 
    }else{
        std::cout << "Invalid integration limits\n";
    }
}

// Simple MC Integration Over Polytopes
Point origin(0);
template 
<
    typename WalkType=BallWalk,
    //typename VolumeType, TODO:: To remove switch interface in volume computation type i.e. for CB / CG / SOB 
    typename Functor,
    typename Polytope
>
void simple_mc_polytope_integrate(Functor Fx,Polytope &P, Uint N,volumeType vType=CB,Point Origin=origin){

    // Polytope volumetric calculation params
    Uint dim = P.dimension();
    int walk_length = 1; 
    NT e=0.1;
    NT volume=0;

    // Check if origin is shifted
    if(Origin.dimension() == 0 ){
        Origin.set_dimension(dim); 
        Origin.set_to_origin();
    }
    
    // Volume calculation for HPolytope
    switch(vType){
        case CB:     
            volume = volume_cooling_balls<BallWalk, RandomNumberGenerator, HPOLYTOPE>(P, e, walk_length).second;
            break;
        case CG: 
            volume = volume_cooling_gaussians<GaussianBallWalk, RandomNumberGenerator>(P, e, walk_length);
            break;
        case SOB: 
            volume = volume_sequence_of_balls<BallWalk, RandomNumberGenerator>(P, e, walk_length);
            break;    
        default:  
            volume = 0;
            std::cout << "Error in enum:vType\n";
            break;
    }    

    std::cout << volume << std::endl;

    // For implementing Uniform Ball Walk Billiards
    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();
    RandomNumberGenerator rng(1);
    Point x0 = inner_ball.first;
    typename WalkType::template Walk<Polytope, RandomNumberGenerator> walk(P,x0,rng);

    NT sum=0;
    for (int i = 0; i < N; i++ ){
        walk.apply(P,x0,walk_length,rng);
        //x0.print();
        sum = sum + Fx(x0+Origin);
    }

    // Final step for integration
    std::cout << "Integral Value over H-Polytope: " << volume * sum / N << "\n";   

}

#endif
