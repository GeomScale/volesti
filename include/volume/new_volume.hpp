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


/////////////////// Random numbers generator
///

template <typename RNGType, typename NT, int ... Ts>
struct BoostRandomNumberGenerator;

template <typename RNGType, typename NT>
struct BoostRandomNumberGenerator<RNGType, NT>
{
    BoostRandomNumberGenerator(int d)
        :   _rng(std::chrono::system_clock::now().time_since_epoch().count())
        ,   _urdist(0, 1)
        ,   _uidist(0, d-1)
        ,   _ndist(0, 1)
    {}

    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<NT> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
    boost::random::normal_distribution<NT> _ndist;
};


template <typename RNGType, typename NT, int Seed>
struct BoostRandomNumberGenerator<RNGType, NT, Seed>
{
    BoostRandomNumberGenerator(int d)
        :   _rng(Seed)
        ,   _urdist(0, 1)
        ,   _uidist(0, d-1)
        ,   _ndist(0, 1)
    {}

    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<NT> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
    boost::random::normal_distribution<NT> _ndist;
};




/////////////////// Random walk helpers
///


template <typename Point>
struct GetDirection
{
    typedef typename Point::FT NT;

    template <typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              RandomNumberGenerator &rng)
    {
        NT normal = NT(0);
        Point p(dim);
        NT* data = p.pointerToData();

        for (unsigned int i=0; i<dim; ++i)
        {
            *data = rng.sample_ndist();
            normal += *data * *data;
            data++;
        }

        normal = NT(1)/std::sqrt(normal);
        p *= normal;

        return p;
    }
};

template <typename Point>
struct GetPointInDsphere
{
    template <typename NT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              NT const& radius,
                              RandomNumberGenerator &rng)
    {
        Point p = GetDirection<Point>::apply(dim, rng);
        NT U = rng.sample_urdist();
        U = std::pow(U, NT(1)/(NT(dim)));
        p *= radius * U;
        return p;
    }
};


/////////////////// Random Walks

// ball walk with uniform target distribution
struct BallWalk
{

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    Walk (Polytope const&, Point&, RandomNumberGenerator&)
    {}

    Walk (BallPolytope const&, Point &, RandomNumberGenerator &)
    {}

    template<typename BallPolytope>
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        const NT delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());

        for (auto j=0u; j<walk_length; ++j)
        {
            Point y = GetPointInDsphere<Point>::apply(P.dimension(),
                                                      delta,
                                                      rng);
            y += p;
            if (P.is_in(y) == -1) p = y;
        }
    }
};

};

// random directions hit-and-run walk with uniform target distribution
struct RDHRWalk
{

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            Point v = GetDirection<Point>::apply(p.dimension(), rng);
            std::pair<NT, NT> bpair = P.line_intersect(_p, v, _lamdas, _Av,
                                                       _lambda);
            _lambda = rng.sample_urdist() * (bpair.first - bpair.second)
                    + bpair.second;
            _p += (_lambda * v);
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           Point &p,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());

        Point v = GetDirection<Point>::apply(p.dimension(), rng);
        std::pair<NT, NT> bpair = P.line_intersect(p, v, _lamdas, _Av);
        _lambda = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
        _p = (_lambda * v) + p;
    }

    Point _p;
    NT _lambda;
    typename Point::Coeff _lamdas;
    typename Point::Coeff _Av;
};

};

// random directions hit-and-run walk with uniform target distribution
struct CDHRWalk
{

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            auto rand_coord_prev = _rand_coord;
            _rand_coord = rng.sample_uidist();
            NT kapa = rng.sample_urdist();
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
    inline void initialize(BallPolytope const& P,
                           Point &p,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _rand_coord = rng.sample_uidist();
        NT kapa = rng.sample_urdist();
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
};

};

// billiard walk for uniform distribution
struct BilliardWalk
{

template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point & p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT diameter = P.ComputeDiameter();
        NT T = rng.sample_urdist() * diameter;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * diameter;
            _v = GetDirection<Point>::apply(n, rng);
            Point p0 = _p;
            int it = 0;
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
    inline void initialize(GenericPolytope const& P,
                           Point &p,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        NT diameter = P.ComputeDiameter();

        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _p = p;
        _v = GetDirection<Point>::apply(n, rng);

        NT T = rng.sample_urdist() * diameter;
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
};

};

///
///////////////////// Random generators' policies

struct PushBackWalkPolicy
{
    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p) const
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

    unsigned int get_nump_PBSmall() const
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
    typename Walk
>
struct RandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator
    >
    static void apply(Polytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng)
    {
        Walk walk(P, p, rng);
        for (unsigned int i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length, rng);
            policy.apply(randPoints, p);
        }
    }
};



////////////////////////////// Algorithms


// ----- VOLUME ------ //

template
<
    typename WalkTypePolicy = CDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, double>,
    typename Polytope
>
double volume_sequence_of_balls(Polytope const& Pin,
                                double const& error = 1.0,
                                unsigned int const& walk_length = 1,
                                unsigned int const& n_threads = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> Ball;
    typedef BallIntersectPolytope<Polytope,Ball> BallPoly;

    typedef typename WalkTypePolicy::template Walk
                                              <
                                                Polytope,
                                                RandomNumberGenerator
                                              > WalkType;
    typedef RandomPointGenerator<WalkType> RandomPointGenerator;

    auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    unsigned int rnum = std::pow(error, -2) * 400 * n * std::log(n);
    RandomNumberGenerator rng(P.dimension());

    //Get the Chebychev ball (largest inscribed ball) with center and radius
    auto InnerBall = P.InnerBall();
    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    VT c_e = Eigen::Map<VT>(&c.get_coeffs()[0], c.dimension());
    P.shift(c_e);
    c=Point(n);

    rnum = rnum/n_threads;
    NT vol = NT(0);

    // Perform the procedure for a number of threads and then take the average
    for (auto t=0u; t<n_threads; t++)
    {
        // Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
#ifdef VOLESTI_DEBUG
        std::cout<<"\nGenerate the first random point in P"<<std::endl;
#endif
        Point p = GetPointInDsphere<Point>::apply(P.dimension(), radius, rng);
        std::list<Point> randPoints; //ds for storing rand points

        PushBackWalkPolicy push_back_policy;
        RandomPointGenerator::apply(P, p, 1, 50*n, randPoints, push_back_policy, rng);

#ifdef VOLESTI_DEBUG
        double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
#endif
        RandomPointGenerator::apply(P, p, rnum-1, walk_length, randPoints,
                                    push_back_policy, rng);

#ifdef VOLESTI_DEBUG
        double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "First random points construction time = "
                  << tstop2 - tstart2 << std::endl;
#endif

        // Construct the sequence of balls
        // a. compute the radius of the largest ball
        NT current_dist, max_dist=NT(0);
        for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit)
        {
            current_dist = (*pit).squared_length();
            if (current_dist > max_dist)
            {
                max_dist=current_dist;
            }
        }
        max_dist = std::sqrt(max_dist);
#ifdef VOLESTI_DEBUG
        std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist
                <<std::endl;
        std::cout<<"\nConstructing the sequence of balls"<<std::endl;
        std::cout<<"---------"<<std::endl;
#endif

        //
        // b. Number of balls
        int nb1 = n * (std::log(radius)/std::log(2.0));
        int nb2 = std::ceil(n * (std::log(max_dist)/std::log(2.0)));

        std::vector<Ball> balls;

        for (auto i=nb1; i<=nb2; ++i)
        {
            if (i == nb1)
            {
                balls.push_back(Ball(c,radius*radius));
                vol = (std::pow(M_PI,n/2.0)*(std::pow(balls[0].radius(), n) ) )
                       / (tgamma(n/2.0+1));
            } else {
                balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
            }
        }
        assert(!balls.empty());

        // Estimate Vol(P)
        typename std::vector<Ball>::iterator bit2=balls.end();
        bit2--;

        while (bit2!=balls.begin())
        {
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
                                        counting_policy, rng);

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
#ifdef VOLESTI_DEBUG
    std::cout<<"rand points = "<<rnum<<std::endl;
    std::cout<<"walk len = "<<walk_length<<std::endl;
    std::cout<<"volume computed: "<<vol<<std::endl;
#endif

    P.free_them_all();
    return vol;
}

#endif
