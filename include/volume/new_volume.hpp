// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NEW_VOLUME_H
#define NEW_VOLUME_H

#include <iterator>
#include <vector>
#include <list>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "ball.h"
#include "ballintersectconvex.h"
#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
#include "gaussian_samplers.h"
#include "gaussian_annealing.h"

#include "khach.h"

/////////////////// Random walk helpers
///
/*
// Pick a random direction as a normilized vector
template <typename RNGType, typename Point, typename NT>
Point get_direction(const unsigned int dim) {

    boost::normal_distribution<> rdist(0,1);
    NT normal = NT(0);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    Point p(dim);
    NT* data = p.pointerToData();

    //RNGType rng2 = var.rng;
    for (unsigned int i=0; i<dim; ++i) {
        *data = rdist(rng);
        normal += *data * *data;
        data++;
    }

    normal=1.0/std::sqrt(normal);
    p *= normal;

    return p;
}


// Pick a random point from a d-sphere
template <typename RNGType, typename Point, typename NT>
Point get_point_on_Dsphere(const unsigned int dim, const NT &radius){
    Point p = get_direction<RNGType, Point, NT>(dim);
    if (radius != 0) p *= radius;
    return p;
}
*/

template <typename RNGType, typename Point>
struct GetDirection
{
    typedef typename Point::FT NT;

    GetDirection(unsigned seed = 1)
        :   _rng(seed)
        ,   _ndist(0,1)
    {}

    Point apply(unsigned int const& dim)
    {
        NT normal = NT(0);
        Point p(dim);
        NT* data = p.pointerToData();

        for (unsigned int i=0; i<dim; ++i)
        {
            *data = _ndist(_rng);
            normal += *data * *data;
            data++;
        }

        normal = NT(1)/std::sqrt(normal);
        p *= normal;

        return p;
    }

private :
    RNGType _rng;
    boost::normal_distribution<> _ndist;
};

template <typename RNGType, typename Point>
struct GetPointInDsphere
{
    GetPointInDsphere(unsigned seed = 1)
        :   _rng(seed)
        ,   _urdist(0,1)
    {}

    template <typename NT>
    Point apply(unsigned int const& dim,
                NT const& radius)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        GetDirection<RNGType, Point> generator(seed);
        Point p = generator.apply(dim);
        NT U = _urdist(_rng);
        U = std::pow(U, NT(1)/(NT(dim)));
        p *= radius * U;
        return p;
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<> _urdist;
};


/////////////////// Random Walks

// ball walk with uniform target distribution
template
<
    typename Polytope,
    typename RNGType
>
struct BallWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    BallWalk(Polytope P)
    {
        P.ComputeInnerBall();
    }

    BallWalk(BallPolytope) {}

    template
    <
        typename BallPolytope
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        const NT delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());

        for (auto j=0; j<walk_length; ++j)
        {
            GetPointInDsphere<RNGType, Point> generator;
            Point y = generator.apply(p.dimension(), delta);
            y += p;
            if (P.is_in(y) == -1) p = y;
        }
    }
};

// random directions hit-and-run walk with uniform target distribution
template
<
    typename Polytope,
    typename RNGType
>
struct RDHRWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    RDHRWalk(Polytope P, unsigned seed = 1)
        :   _rng(seed)
        ,   _urdist(0,1)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    RDHRWalk(BallPolytope P, unsigned seed = 1)
        :   _rng(seed)
        ,   _urdist(0,1)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    template
    <
        typename BallPolytope
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        for (auto j=0; j<walk_length; ++j)
        {
            Point v = get_direction<RNGType, Point, NT>(p.dimension());
            std::pair<NT, NT> bpair = P.line_intersect(_p, v, _lamdas, _Av,
                                                       _lambda);
            _lambda = _urdist(_rng) * (bpair.first - bpair.second)
                    + bpair.second;
            _p += (_lambda * v);
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    void initialize(BallPolytope &P, Point &p)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());

        Point v = get_direction<RNGType, Point, NT>(p.dimension());
        std::pair<NT, NT> bpair = P.line_intersect(p, v, _lamdas, _Av);
        _lambda = _urdist(_rng) * (bpair.first - bpair.second) + bpair.second;
        _p = (_lambda * v) + p;
    }

    Point _p;
    NT _lambda;
    typename Point::Coeff _lamdas;
    typename Point::Coeff _Av;
    RNGType _rng;
    boost::random::uniform_real_distribution<> _urdist;
};


// random directions hit-and-run walk with uniform target distribution
template
<
    typename Polytope,
    typename RNGType
>
struct CDHRWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    CDHRWalk(Polytope P, unsigned seed = 1)
        :   _rng(seed)
        ,   _urdist(0,1)
        ,   _uidist(0, P.dimension()-1)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    CDHRWalk(BallPolytope P, unsigned seed = 1)
        :   _rng(seed)
        ,   _urdist(0,1)
        ,   _uidist(0, P.dimension()-1)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    template
    <
        typename BallPolytope
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        for (auto j=0; j<walk_length; ++j)
        {
            auto rand_coord_prev = _rand_coord;
            _rand_coord = _uidist(_rng);
            NT kapa = _urdist(_rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(_p,
                                                             _p_prev,
                                                             _rand_coord,
                                                             rand_coord_prev,
                                                             _lamdas);
            _p_prev = _p;
            _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                         * (bpair.second - bpair.first));
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    void initialize(BallPolytope &P, Point &p)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _rand_coord = _uidist(_rng);
        NT kapa = _urdist(_rng);
        _p = p;
        std::pair<NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord, _lamdas);
        _p_prev = _p;
        _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                    * (bpair.second - bpair.first));
    }

    unsigned int _rand_coord;
    Point _p;
    Point _p_prev;
    typename Point::Coeff _lamdas;
    RNGType _rng;
    boost::random::uniform_real_distribution<> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
};


// billiard walk for uniform distribution
template
<
    typename Polytope,
    typename RNGType
>
struct BilliardWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    BilliardWalk(Polytope P, unsigned seed = 1)
        :   _rng(seed)
        ,   _urdist(0,1)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    BilliardWalk(BallPolytope P, unsigned seed = 1)
        :   _rng(seed)
        ,   _urdist(0,1)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    template
    <
        typename GenericPolytope
    >
    void apply(GenericPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        unsigned int n = P.dimension();
        NT diameter = P.ComputeDiameter();
        NT T = _urdist(_rng) * diameter;
        const NT dl = 0.995;

        for (auto j=0; j<walk_length; ++j)
        {
            T = _urdist(_rng) * diameter;
            _v = get_direction<RNGType, Point, NT>(n);
            Point p0 = _p;
            auto it = 0;
            while (it < 10*n)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == 10*n) _p = p0;
        }
        p = _p;
    }

private :

    template
    <
        typename GenericPolytope
    >
    void initialize(GenericPolytope &P,
                    Point &p)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        NT diameter = P.ComputeDiameter();

        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _p = p;
        _v = get_direction<RNGType, Point, NT>(n);

        NT T = _urdist(_rng) * diameter;
        std::pair<NT, int> pbpair
                = P.line_positive_intersect(_p, _v, _lambdas, _Av);
        if (T <= pbpair.first) {
            _p += (T * _v);
            _lambda_prev = T;
            return;
        }
        _lambda_prev = dl * pbpair.first;
        _p += (_lambda_prev * _v);
        T -= _lambda_prev;
        P.compute_reflection(_v, _p, pbpair.second);
    }

    Point _p;
    Point _v;
    NT _lambda_prev;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
    RNGType _rng;
    boost::random::uniform_real_distribution<> _urdist;
};


///
/// Random generators' policies

struct PushBackWalkPolicy
{
    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p)
    {
        randPoints.push_back(p);
    }
};

template <typename BallPoly>
struct CountingWalkPolicy
{
    CountingWalkPolicy(unsigned int const& nump_PBSmall, BallPoly const& PBSmall)
        :   _nump_PBSmall(nump_PBSmall)
        ,   _PBSmall(PBSmall)
    {}

    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p)
    {
        if (_PBSmall.second().is_in(p) == -1)//is in
        {
            randPoints.push_back(p);
            ++_nump_PBSmall;
        }
    }

    unsigned int get_nump_PBSmall()
    {
        return _nump_PBSmall;
    }

private :
    unsigned int _nump_PBSmall;
    BallPoly _PBSmall;
};


////////////////////////////// Random Point Generators
///

template
<
    typename Walk,
    typename RNGType
>
struct RandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename PointList,
        typename WalkPolicy
    >
    static void apply(Polytope &P,
                      Point &p,   // a point to start
                      const unsigned int rnum,
                      const unsigned int walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        Walk walk(P);
        for (auto i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length);
            policy.apply(randPoints, p);
        }
    }
};



////////////////////////////// Algorithms


// ----- VOLUME ------ //

template
<
    typename Polytope,
    typename RNGType = boost::mt19937,
    typename WalkType = BallWalk<Polytope,RNGType>
>
double volume(Polytope &P,
              double error = 1.0,
              unsigned int walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> Ball;
    typedef BallIntersectPolytope<Polytope,Ball> BallPoly;

    typedef RandomPointGenerator<WalkType, RNGType> RandomPointGenerator;

    unsigned int n = P.dimension();
    unsigned int rnum = std::pow(error, -2) * 400 * n * std::log(n);
    unsigned int n_threads = 1;

    //Get the Chebychev ball (largest inscribed ball) with center and radius
    auto InnerBall = P.ComputeInnerBall();
    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    VT c_e = Eigen::Map<VT>(&c.get_coeffs()[0], c.dimension());
    P.shift(c_e);
    c=Point(n);

    rnum=rnum/n_threads;
    NT vol = NT(0);

    // Perform the procedure for a number of threads and then take the average
    for(unsigned int t=0; t<n_threads; t++){
        // Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
#ifdef VOLESTI_DEBUG
        std::cout<<"\nGenerate the first random point in P"<<std::endl;
#endif
        Point p = get_point_on_Dsphere<RNGType , Point>(n, radius);
        std::list<Point> randPoints; //ds for storing rand points

        PushBackWalkPolicy push_back_policy;
        RandomPointGenerator::apply(P, p, 1, 50*n, randPoints, push_back_policy);

        // Sample "rnum" points from P
#ifdef VOLESTI_DEBUG
        double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
#endif
        RandomPointGenerator::apply(P, p, rnum-1, walk_length, randPoints,
                                    push_back_policy);

#ifdef VOLESTI_DEBUG
        double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "First random points construction time = "
                  << tstop2 - tstart2 << std::endl;
#endif

        // Construct the sequence of balls
        // a. compute the radius of the largest ball
        NT current_dist, max_dist=NT(0);
        for(typename  std::list<Point>::iterator pit=randPoints.begin();
            pit!=randPoints.end(); ++pit){
            current_dist=(*pit).squared_length();
            if(current_dist>max_dist){
                max_dist=current_dist;
            }
        }
        max_dist=std::sqrt(max_dist);
#ifdef VOLESTI_DEBUG
        std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist
                <<std::endl;
        std::cout<<"\nConstructing the sequence of balls"<<std::endl;
#endif

        //
        // b. Number of balls
        int nb1 = n * (std::log(radius)/std::log(2.0));
        int nb2 = std::ceil(n * (std::log(max_dist)/std::log(2.0)));

        std::vector<Ball> balls;

        for (int i=nb1; i<=nb2; ++i)
        {
            if (i==nb1)
            {
                balls.push_back(Ball(c,radius*radius));
                vol = (std::pow(M_PI,n/2.0)*(std::pow(balls[0].radius(), n) ) )
                       / (tgamma(n/2.0+1));
            } else {
                balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
            }
        }
        assert(!balls.empty());

#ifdef VOLESTI_DEBUG
        std::cout<<"---------"<<std::endl;
#endif

        // Estimate Vol(P)

        typename std::vector<Ball>::iterator bit2=balls.end();
        bit2--;

        while(bit2!=balls.begin()){

            //each step starts with some random points in PBLarge stored
            //in list "randPoints", these points have been generated in a
            //previous step

            BallPoly PBLarge(P,*bit2);
            --bit2;
            BallPoly PBSmall(P,*bit2);

#ifdef VOLESTI_DEBUG
            std::cout<<"("<<balls.end()-bit2<<"/"<<balls.end()-balls.begin()
                     <<")Ball ratio radius="
                     <<PBLarge.second().radius()<<","
                     <<PBSmall.second().radius()<<std::endl;
            std::cout<<"Points in PBLarge="<<randPoints.size()<<std::endl;
#endif

            // choose a point in PBLarge to be used to generate more rand points
            Point p_gen = *randPoints.begin();

            // num of points in PBSmall and PBLarge
            unsigned int nump_PBSmall = 0;
            unsigned int nump_PBLarge = randPoints.size();

            //keep the points in randPoints that fall in PBSmall
            typename std::list<Point>::iterator rpit=randPoints.begin();
            while (rpit!=randPoints.end())
            {
                if (PBSmall.second().is_in(*rpit) == 0)//not in
                {
                    rpit=randPoints.erase(rpit);
                } else {
                    ++nump_PBSmall;
                    ++rpit;
                }
            }

#ifdef VOLESTI_DEBUG
            std::cout<<"Points in PBSmall="<<randPoints.size()
                     <<"\nRatio= "<<NT(nump_PBLarge)/NT(nump_PBSmall)
                     <<std::endl;
            std::cout<<"Generate "<<rnum-nump_PBLarge<<" more "<<std::endl;
#endif

            CountingWalkPolicy<BallPoly> counting_policy(nump_PBSmall, PBSmall);
            RandomPointGenerator::apply(PBLarge, p_gen, rnum-nump_PBLarge,
                             walk_length, randPoints,
                             counting_policy);

            nump_PBSmall = counting_policy.get_nump_PBSmall();

            vol *= NT(rnum)/NT(nump_PBSmall);

#ifdef VOLESTI_DEBUG
            std::cout<<nump_PBSmall<<"/"<<rnum<<" = "<<NT(rnum)/nump_PBSmall
                     <<"\ncurrent_vol = "<<vol
                     <<"\n--------------------------"<<std::endl;
#endif

            //don't continue in pairs of balls that are almost inside P, i.e. ratio ~= 2
        }
    }
    //vol=round_value*vol;
#ifdef VOLESTI_DEBUG
    std::cout<<"rand points = "<<rnum<<std::endl;
    std::cout<<"walk len = "<<walk_length<<std::endl;
    //std::cout<<"round_value: "<<round_value<<std::endl;
    std::cout<<"volume computed: "<<vol<<std::endl;
#endif

    P.free_them_all();
    return vol;
}

#endif
