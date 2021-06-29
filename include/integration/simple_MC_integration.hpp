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
typedef std::vector<Point> Points;
typedef HPolytope<Point> HPOLYTOPE;
typedef boost::mt19937 RNGType;
typedef BoostRandomNumberGenerator<RNGType, NT> RandomNumberGenerator;
typedef typename HPolytope<Point>::MT MT;
typedef typename HPolytope<Point>::VT VT; 

typedef const unsigned int Uint;  // positive constant value for no of samples & dimensions
enum volType { CB , CG , SOB }; // Volume type for polytope
const Point pt(0); // Required to initialize points in function params
typedef std::pair<Point, Point> Limits;

// To check if two n-dimensional points ensure valid limits in integration 
bool validLimits(Point LL, Point UL){
    if( UL.dimension() == LL.dimension() ) {
        for ( int i = 0 ; i< LL.dimension() ; i++){
            if( LL[i] > UL[i]){
            std::cerr << "Invalid integration limits\n";
            return false;
            }
        }
        return true;
    }else{
        std::cerr << "Invalid integration limits\n";
        return false;
    }
}

// Initialize Limit Point
Point initLimit(NT L[], Uint dim){
    Point pt(dim);
    for(int i=0 ; i<dim ; i++){
        pt.set_coord(i, L[i]);
    }
    return pt;
}

//Initialize to [-1,1]^n
void initiateLimits(Point& LL, Point& UL, int dim){
    LL.set_dimension(dim);
    UL.set_dimension(dim);
    LL.set_to_origin();
    UL.set_to_origin();
    for(int i=0 ; i<dim ; i++){
        LL.set_coord(i,-1);
        UL.set_coord(i,1);
    }
}

// Hyper-Rectangle volume calculator in n-dimensions
NT hyper_rect_volume(Point LL, Point UL){
    NT product=1;
    if(validLimits(LL,UL)){
        for(int i=0; i<LL.dimension(); ++i){
            product = product * ( UL[i] - LL[i] );
        }
        return product;
    }
    else return 0;
}

// To sample a point between two n-dimensional points using inbuilt random sampling
Point samplerBWLimits(Point LL, Point UL){
    Point sample_point(LL.dimension());
    for(int i=0; i<LL.dimension(); ++i){
        sample_point.set_coord(i , LL[i] + (NT)(rand()) / ((NT)(RAND_MAX/(UL[i] - LL[i]))) );
    }
    return sample_point;
}

// Simple MC Integration over Hyper-Rectangles
template
<
    typename WalkType=BallWalk,
    typename Functor
>
NT simple_mc_integrate(Functor Fx, Uint dim, Uint N ,Point LL=pt ,Point UL=pt ,int walk_length=10, NT e=0.1){
 
    NT sum = 0;
    if(validLimits(LL,UL)){

        if(LL.dimension()==0 && UL.dimension()==0 ) initiateLimits(LL, UL, dim);

        // Creating an MT & VT for HPolytope(Hyper-Rectangle) for integration limits using LL & UL
        MT mt(dim*2,dim); 
        mt.setZero();
        VT vt(dim*2); 
        vt.setZero();

        for(int i=0 ; i<dim ; i++){
            mt(i,i)=1;
            vt(i)= UL[i];
            mt(dim+i,i)=-1;
            vt(dim+i)=LL[i]*(-1);
        }

        // Initialization of H-Polytope
        HPOLYTOPE P(dim,mt,vt);
        // P.print();
    
        // Setting up params for random walks
        std::pair<Point, NT> inner_ball = P.ComputeInnerBall();
        RandomNumberGenerator rng(1);
        Point x0 = inner_ball.first;
        typename WalkType::template Walk<HPOLYTOPE, RandomNumberGenerator> walk(P, x0, rng);
        
        for (int i = 0; i <=N; i++ ) {
            walk.apply(P, x0, walk_length, rng);
            sum += Fx(x0);
        }    

        NT volume = hyper_rect_volume(LL,UL);
        // std::cout << volume << std::endl;

        NT integration_value = volume * sum / N ;
        return integration_value;
    }else{
        std::cerr << "Invalid integration limits\n";
        return 0;
    }
}


// Simple MC Integration Over Polytopes
template 
<
    typename WalkType=BallWalk,
    typename Polytope=HPOLYTOPE,
    typename RNG=RandomNumberGenerator,
    typename Functor
>
NT simple_mc_polytope_integrate(Functor Fx, Uint dim, Polytope &P, Uint N=10000, volType vT=SOB, int walk_length=1, NT e=0.1 ,Point Origin=pt){

    // P.print();

    // Check if origin is shifted
    if(Origin.dimension() == 0 ){
        Origin.set_dimension(dim); 
        Origin.set_to_origin();
    }else if(Origin.dimension() != dim){
        std::cerr << "Enter Polytope sample point shift properly" << std::endl;
        return 0;
    }

    // Volume calculation for HPolytope
    NT volume=0;
    
    switch(vT){
    case CB:     
        volume = volume_cooling_balls<BallWalk, RNG, Polytope>(P, e, walk_length).second; 
        break;
    case CG: 
        volume = volume_cooling_gaussians<GaussianBallWalk, RNG, Polytope>(P, e, walk_length);
        break;
    case SOB: 
        volume = volume_sequence_of_balls<BallWalk, RNG, Polytope>(P, e, walk_length);
        break;
    default:
        std::cerr << "Error in volume type: CB / SOB / CG" << std::endl;
        volume = 0;
        break;
    }

    // Volume of the Polytope
    // std::cout << "Volume of the Polytope = " << volume << std::endl;

    // For implementing Uniform Walks
    RNG rng(1);
    std::pair<Point, NT> inner_ball = P.ComputeInnerBall();
    Point x0 = inner_ball.first;
    typename WalkType::template Walk<Polytope,RNG> walk(P, x0, rng);

    // For storing sampled points
    Points points; 
    NT sum=0;

    // Applying and walking through Uniform Walks + Storing Points in Vector<Point>
    for (int i = 0; i < N; i++ ){
        walk.apply(P, x0, walk_length, rng);
        sum+=Fx(x0+Origin);

        // points.push_back(x0+Origin);
        // (x0+Origin).print();
    }

    // To print sampled points
    // for(auto i_th_point : points){
    //     i_th_point.print();
    // }

    /*
    Core idea of Monte Carlo Integration algorithm used here : https://en.wikipedia.org/wiki/Monte_Carlo_integration#Overview
    */

    // Integration Value
    NT integration_value = volume * sum / N ;
    return integration_value;
}

#endif
