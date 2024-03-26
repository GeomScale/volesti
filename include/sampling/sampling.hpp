// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.


#ifndef SAMPLE_ONLY_H
#define SAMPLE_ONLY_H
template <typename WalkTypePolicy,
          typename PointList,
          typename Polytope,
          typename RandomNumberGenerator,
          typename Point>
void uniform_sampling(PointList &randPoints,
                      Polytope &P,
                      RandomNumberGenerator &rng,
                      const unsigned int &walk_len,
                      const unsigned int &rnum,
                      Point &starting_point,  
                      unsigned int const& nburns) {

    typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk(P, starting_point, rng);

    Point p = starting_point; 

    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walk.apply(P, p, walk_len, rng);  
        }
    }

    for (unsigned int i = 0; i < rnum; ++i) {
        walk.apply(P, p, walk_len, rng); 
        randPoints.push_back(p);  
    }
}



template <
        typename PointList,
        typename Polytope,
        typename RandomNumberGenerator,
        typename WalkTypePolicy,
        typename Point
>
void uniform_sampling(PointList &randPoints,
                      Polytope &P,
                      RandomNumberGenerator &rng,
                      WalkTypePolicy &WalkType, 
                      const unsigned int &walk_len,
                      const unsigned int &rnum,
                      Point &starting_point,  
                      unsigned int const& nburns) {

    typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk(P, starting_point, rng);

    Point p = starting_point; 
    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walk.apply(P, p, walk_len, rng);  
        }
    }

    randPoints.clear();

    for (unsigned int i = 0; i < rnum; ++i) {
        walk.apply(P, p, walk_len, rng);  
        randPoints.push_back(p);  
    }
}



template <
    typename WalkTypePolicy,
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename Point
>
void uniform_sampling_boundary(PointList &randPoints,
                               Polytope &P,
                               RandomNumberGenerator &rng,
                               const unsigned int &walk_len,
                               const unsigned int &rnum,
                               Point &starting_point, 
                               unsigned int const& nburns) {

    
    typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk(P, starting_point, rng);

    Point p1 = starting_point; 

    Point p2 = starting_point; 
    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i){
            walk.apply(P, p1, p2, walk_len, rng);
        }
    }

    
    randPoints.clear();

    for (unsigned int i = 0; i < rnum; ++i) {
        walk.apply(P, p1, p2, walk_len, rng); 
        randPoints.push_back(p1); 
    }
}


template <
    typename WalkTypePolicy,
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename NT,
    typename Point
>
void gaussian_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const NT &a,
                       const Point &starting_point, 
                       unsigned int const& nburns) {

    typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk(P, starting_point, a, rng);

    Point p = starting_point; 
    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walk.apply(P, p, a, walk_len, rng); 
        }
    }

    randPoints.clear(); 

    for (unsigned int i = 0; i < rnum; ++i) {
        walk.apply(P, p, a, walk_len, rng); 
        randPoints.push_back(p); 
    }
}


template <
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename WalkTypePolicy,
    typename NT,
    typename Point
>
void gaussian_sampling(PointList &randPoints,
                       Polytope &P,
                       RandomNumberGenerator &rng,
                       WalkTypePolicy &WalkType, 
                       const unsigned int &walk_len,
                       const unsigned int &rnum,
                       const NT &a,
                       const Point &starting_point,
                       unsigned int const& nburns) {

    
    typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk(P, starting_point, a, rng);

    Point p = starting_point; 

    
    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walk.apply(P, p, a, walk_len, rng); 
        }
    }

    
    randPoints.clear();

    
    for (unsigned int i = 0; i < rnum; ++i) {
        walk.apply(P, p, a, walk_len, rng); 
        randPoints.push_back(p); 
    }
}


template <
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename WalkTypePolicy,
    typename NT,
    typename Point,
    typename NegativeGradientFunctor,
    typename NegativeLogprobFunctor,
    typename Solver
>
void logconcave_sampling(PointList &randPoints,
                         Polytope &P,
                         RandomNumberGenerator &rng,
                         const unsigned int &walk_len,
                         const unsigned int &rnum,
                         const Point &starting_point,
                         unsigned int const& nburns,
                         NegativeGradientFunctor &F,
                         NegativeLogprobFunctor &f) {
    
    typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver> walk(P, starting_point, F, f, rng);

    Point p = starting_point; 

    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walk.apply(P, p, walk_len, rng); 
        }
    }

    randPoints.clear();

    for (unsigned int i = 0; i < rnum; ++i) {
        walk.apply(P, p, walk_len, rng); 
        randPoints.push_back(p); 
    }
}

#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
template <
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename WalkTypePolicy,
    typename NT,
    typename Point,
    typename NegativeGradientFunctor,
    typename NegativeLogprobFunctor,
    typename HessianFunctor,
    typename Solver
>
void crhmc_sampling(PointList &randPoints,
                    Polytope &P,
                    RandomNumberGenerator &rng,
                    const int walk_len,
                    const unsigned int rnum,
                    const unsigned int nburns,
                    NegativeGradientFunctor &F,
                    NegativeLogprobFunctor &f,
                    HessianFunctor &h,
                    int simdLen = 1,
                    bool raw_output=false) {

    typedef typename Polytope::MT MatrixType;
    typedef crhmc_input<MatrixType, Point, NegativeLogprobFunctor, NegativeGradientFunctor, HessianFunctor> Input;

    Input input = convert2crhmc_input<Input, Polytope, NegativeLogprobFunctor, NegativeGradientFunctor, HessianFunctor>(P, f, F, h);
    typedef crhmc_problem<Point, Input> CrhmcProblem;
    CrhmcProblem problem = CrhmcProblem(input);

    if(problem.terminate) { return; }

    typename WalkTypePolicy::template parameters<NT, NegativeGradientFunctor> crhmc_params(F, P.dimension(), problem.options);
    Point p = Point(problem.center);

    typename WalkTypePolicy::template Walk<Point, CrhmcProblem, RandomNumberGenerator, NegativeGradientFunctor, NegativeLogprobFunctor, Solver>
    crhmc_walk(problem, p, F, f, crhmc_params);

    for(unsigned int i = 0; i < nburns; ++i) {
        crhmc_walk.apply(rng, walk_len);
    }

    for(unsigned int i = 0; i < rnum; ++i) {
        crhmc_walk.apply(rng, walk_len);
        randPoints.push_back(crhmc_walk.getPoint()); 
    }
}







#include "ode_solvers/ode_solvers.hpp"
template <
        typename Polytope,
        typename RNGType,
        typename PointList,
        typename NegativeGradientFunctor,
        typename NegativeLogprobFunctor,
        typename HessianFunctor,
        typename CRHMCWalk,
        int simdLen=1
>
void execute_crhmc(Polytope &P, RNGType &rng, PointList &randPoints,
                  unsigned int const& walkL, unsigned int const& numpoints,
                  unsigned int const& nburns, NegativeGradientFunctor *F=NULL,
                  NegativeLogprobFunctor *f=NULL, HessianFunctor *h=NULL, bool raw_output= false){
typedef typename Polytope::MT MatrixType;
typedef typename Polytope::PointType Point;
typedef typename Point::FT NT;
if(h!=NULL){
typedef  crhmc_input
  <
    MatrixType,
    Point,
    NegativeLogprobFunctor,
    NegativeGradientFunctor,
    HessianFunctor
  > Input;
typedef crhmc_problem<Point, Input> CrhmcProblem;
crhmc_sampling <
  PointList,
  Polytope,
  RNGType,
  CRHMCWalk,
  NT,
  Point,
  NegativeGradientFunctor,
  NegativeLogprobFunctor,
  HessianFunctor,
  ImplicitMidpointODESolver <
  Point,
  NT,
  CrhmcProblem,
  NegativeGradientFunctor,
  simdLen
  >
>(randPoints, P, rng, walkL, numpoints, nburns, *F, *f, *h, simdLen, raw_output);
}else{
  typedef  crhmc_input
        <
                MatrixType,
                Point,
                NegativeLogprobFunctor,
                NegativeGradientFunctor,
                ZeroFunctor<Point>
        > Input;
  typedef crhmc_problem<Point, Input> CrhmcProblem;
  ZeroFunctor<Point> zerof;
crhmc_sampling <
  PointList,
  Polytope,
  RNGType,
  CRHMCWalk,
  NT,
  Point,
  NegativeGradientFunctor,
  NegativeLogprobFunctor,
  ZeroFunctor<Point>,
  ImplicitMidpointODESolver <
  Point,
  NT,
  CrhmcProblem,
  NegativeGradientFunctor,
  simdLen
  >
>(randPoints, P, rng, walkL, numpoints, nburns, *F, *f, zerof, simdLen, raw_output);
}
}




template <
    typename WalkTypePolicy,
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename NT,
    typename Point
>
void exponential_sampling(
    PointList &randPoints,
    Polytope &P,
    RandomNumberGenerator &rng,
    const unsigned int &walk_len,
    const unsigned int &rnum,
    const Point &c,
    const NT &a,
    Point &starting_point, 
    unsigned int const& nburns
) {
    typedef typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk;

    PushBackWalkPolicy push_back_policy;
    Point p = starting_point; 

    typedef ExponentialRandomPointGenerator<walk> RandomPointGenerator;

    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walk tmpWalk(P, p, rng); 
            tmpWalk.apply(p, c, a, walk_len);  
        }
    }

    for (unsigned int i = 0; i < rnum; ++i) {
        walk tmpWalk(P, p, rng);  
        tmpWalk.apply(p, c, a, walk_len);  
        randPoints.push_back(p);  
    }
}



template <
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename WalkTypePolicy,
    typename NT,
    typename Point
>
void exponential_sampling(
    PointList &randPoints,
    Polytope &P,
    RandomNumberGenerator &rng,
    WalkTypePolicy &WalkType, 
    const unsigned int &walk_len,
    const unsigned int &rnum,
    const Point &c,
    const NT &a,
    Point &starting_point, 
    unsigned int const& nburns
) {
    typedef typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk;
    walk walkInstance(P, starting_point, rng, WalkType.param);  

    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walkInstance.apply(starting_point, c, a, walk_len);  
        }
    }
    for (unsigned int i = 0; i < rnum; ++i) {
        walkInstance.apply(starting_point, c, a, walk_len);  
        randPoints.push_back(starting_point);  
    }
}



template <
    typename WalkTypePolicy,
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename NT,
    typename Point
>
void exponential_sampling(
    PointList &randPoints,
    Polytope &P,
    RandomNumberGenerator &rng,
    const unsigned int &walk_len,
    const unsigned int &rnum,
    const Point &c,
    const NT &a,
    const NT &eta,
    Point &starting_point,  
    unsigned int const& nburns
) {
    typedef typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk;

    walk walkInstance(P, starting_point, rng);
    if (nburns > 0) {
        for (unsigned int i = 0; i < nburns; ++i) {
            walkInstance.apply(starting_point, c, a, eta, walk_len);  
        }
    }

    
    randPoints.clear(); 
    for (unsigned int i = 0; i < rnum; ++i) {
        walkInstance.apply(starting_point, c, a, eta, walk_len); 
        randPoints.push_back(starting_point);  
    }
}



template <
    typename PointList,
    typename Polytope,
    typename RandomNumberGenerator,
    typename WalkTypePolicy,
    typename NT,
    typename Point
>
void exponential_sampling(PointList &randPoints,
                          Polytope &P,
                          RandomNumberGenerator &rng,
                          WalkTypePolicy &WalkType,
                          const unsigned int &walk_len,
                          const unsigned int &rnum,
                          const Point &c,
                          const NT &a,
                          const NT &eta, 
                          const Point &starting_point,
                          unsigned int const& nburns)
{
    
    typename WalkTypePolicy::template Walk<Polytope, RandomNumberGenerator> walk(P, starting_point, rng, WalkType.param);

    Point p = starting_point;

    for(unsigned int i = 0; i < nburns; ++i) {
        walk.apply(p, walk_len, c, a, eta, rng); 
    }

    randPoints.clear();

    for(unsigned int i = 0; i < rnum; ++i) {
        walk.apply(p, walk_len, c, a, eta, rng); 
        randPoints.push_back(p); 
    }
}


#endif
