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

/////////////////// Random Walks

// ball walk with uniform target distribution
template <typename Point>
struct BallWalk
{
    typedef typename Point::FT NT;

    template
    <
        typename RNGType,
        typename Polytope,
        typename PointList,
        typename WalkPolicy
    >
    void apply(Polytope &P,
               Point &p,   // a point to start
               const unsigned int rnum,
               const unsigned int walk_length,
               PointList &randPoints,
               WalkPolicy &policy)
    {
        const NT delta = (P.InnerBall()).second;

        for (auto i=0; i<rnum; ++i)
        {
            for (auto j=0; j<walk_length; ++j)
            {
                Point y = get_point_in_Dsphere<RNGType, Point>(p.dimension(), delta);
                y = y + p;
                if (P.is_in(y)==-1) p = y;
            }
            policy.apply(randPoints, p);
        }
    }
};

// random directions hit-and-run walk with uniform target distribution
template <typename Point>
struct RDHRWalk
{
    typedef typename Point::FT NT;

    template
    <
        typename RNGType,
        typename Polytope,
        typename PointList,
        typename WalkPolicy
    >
    void apply(Polytope &P,
               Point &p,   // a point to start
               const unsigned int rnum,
               const unsigned int walk_length,
               PointList &randPoints,
               WalkPolicy &policy)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        initialize<RNGType>(P, p);

        for (auto i=0; i<rnum; ++i)
        {
            for (auto j=0; j<walk_length; ++j)
            {
                Point v = get_direction<RNGType, Point, NT>(p.dimension());
                std::pair<NT, NT> bpair = P.line_intersect(_p, v, _lamdas, _Av,
                                                           _lambda);
                _lambda = urdist(rng) * (bpair.first - bpair.second)
                        + bpair.second;
                _p = (_lambda * v) + _p;
            }
            policy.apply(randPoints, _p);
        }
    }

private :

    template <typename RNGType, typename Polytope>
    void initialize(Polytope &P, Point &p)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        _lamdas.resize(P.num_of_hyperplanes(), NT(0));
        _Av.resize(P.num_of_hyperplanes(), NT(0));
        Point v = get_direction<RNGType, Point, NT>(p.dimension());
        std::pair<NT, NT> bpair = P.line_intersect(p, v, _lamdas, _Av);
        _lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        _p = (_lambda * v) + p;
    }

    Point _p;
    NT _lambda;
    std::vector<NT> _lamdas;
    std::vector<NT> _Av;
};


// random directions hit-and-run walk with uniform target distribution
template <typename Point>
struct CDHRWalk
{
    typedef typename Point::FT NT;

    template
    <
        typename RNGType,
        typename Polytope,
        typename PointList,
        typename WalkPolicy
    >
    void apply(Polytope &P,
               Point &p,   // a point to start
               const unsigned int rnum,
               const unsigned int walk_length,
               PointList &randPoints,
               WalkPolicy &policy)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, p.dimension()-1);

        initialize<RNGType>(P, p);

        for (auto i=0; i<rnum; ++i)
        {
            for (auto j=0; j<walk_length; ++j)
            {
                auto rand_coord_prev = _rand_coord;
                _rand_coord = uidist(rng);
                NT kapa = urdist(rng);
                std::pair<NT, NT> bpair = P.line_intersect_coord(_p,
                                                                 _p_prev,
                                                                 _rand_coord,
                                                                 rand_coord_prev,
                                                                 _lamdas);
                _p_prev = _p;
                _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                             * (bpair.second - bpair.first));
            }
            policy.apply(randPoints, _p);
        }
    }

private :

    template <typename RNGType, typename Polytope>
    void initialize(Polytope &P, Point &p)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, p.dimension()-1);

        _lamdas.resize(P.num_of_hyperplanes(), NT(0));
        _rand_coord = uidist(rng);
        NT kapa = urdist(rng);
        _p=p;
        std::pair<NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord, _lamdas);
        _p_prev = _p;
        _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                    * (bpair.second - bpair.first));
    }

    unsigned int _rand_coord;
    Point _p;
    Point _p_prev;
    std::vector<NT> _lamdas;
};


// billiard walk for uniform distribution
template <typename Point>
struct BilliardWalk
{
    typedef typename Point::FT NT;

    template
    <
        typename RNGType,
        typename Polytope,
        typename PointList,
        typename WalkPolicy
    >
    void apply(Polytope &P,
               Point &p,   // a point to start
               const unsigned int rnum,
               const unsigned int walk_length,
               PointList &randPoints,
               WalkPolicy &policy)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        unsigned int n = P.dimension();
        NT T = urdist(rng) * P.ComputeDiameter();
        const NT dl = 0.995;
        NT diameter = P.ComputeDiameter();

        _lambdas.resize(P.num_of_hyperplanes(), NT(0));
        _Av.resize(P.num_of_hyperplanes(), NT(0));

        _p = p;

        for (auto i=0; i<rnum; ++i)
        {
            for (auto j=0; j<walk_length; ++j)
            {
                T = urdist(rng) * diameter;
                _v = get_direction<RNGType, Point, NT>(n);
                Point p0 = _p;

                auto it = 0;
                while (it < 10*n)
                {
                    std::pair<NT, int> pbpair;
                    if (i==0 && j==0)
                    {
                        pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av);
                    } else {
                        pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev);
                    }
                    if (T <= pbpair.first) {
                        _p = (T * _v) + _p;
                        _lambda_prev = T;
                        break;
                    }

                    _lambda_prev = dl * pbpair.first;
                    _p = (_lambda_prev * _v) + _p;
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, pbpair.second);
                    it++;
                }
                if (it == 10*n) _p = p0;
            }
            policy.apply(randPoints, _p);
        }
    }

private :
    unsigned int _rand_coord;
    Point _p;
    Point _v;
    NT _lambda_prev;
    std::vector<NT> _lambdas;
    std::vector<NT> _Av;
};



////////////////////////////// Algorithms

#include "khach.h"


// ----- ROUNDING ------ //
// main rounding function
template <typename Polytope, typename Point, typename Parameters, typename NT>
std::pair<NT, NT> rounding_min_ellipsoid2(Polytope &P,
                                         const std::pair<Point,NT> &InnerBall,
                                         const Parameters &var)
{
    typedef typename Polytope::MT 	MT;
    typedef typename Polytope::VT 	VT;
    typedef typename Parameters::RNGType RNGType;
    unsigned int n=var.n, walk_len=var.walk_steps, i, j = 0;
    Point c = InnerBall.first;
    NT radius = InnerBall.second;
    std::list<Point> randPoints; //ds for storing rand points
    if (!P.get_points_for_rounding(randPoints)) {  // If P is a V-polytope then it will store its vertices in randPoints
        // If P is not a V-Polytope or number_of_vertices>20*domension
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        Point p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        p = p + c;

        //use a large walk length e.g. 1000
        rand_point_generator(P, p, 1, 10*n, randPoints, var);
        // 3. Sample points from P
        unsigned int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
        randPoints.clear();
        if (var.bill_walk) {
            rand_point_generator(P, p, num_of_samples, 5, randPoints, var);
        } else {
            rand_point_generator(P, p, num_of_samples, 10 + n / 10, randPoints, var);
        }
    }

    // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
    boost::numeric::ublas::matrix<double> Ap(n,randPoints.size());
    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            Ap(i,j)=double(*qit);
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n);
    boost::numeric::ublas::vector<double> c2(n);
    size_t w=1000;
    KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

    MT E(n,n);
    VT e(n);

    //Get ellipsoid matrix and center as Eigen objects
    for(unsigned int i=0; i<n; i++){
        e(i)=NT(c2(i));
        for (unsigned int j=0; j<n; j++){
            E(i,j)=NT(Q(i,j));
        }
    }


    //Find the smallest and the largest axes of the elliposoid
    Eigen::EigenSolver<MT> eigensolver(E);
    NT rel = std::real(eigensolver.eigenvalues()[0]);
    NT Rel = std::real(eigensolver.eigenvalues()[0]);
    for(unsigned int i=1; i<n; i++){
        if(std::real(eigensolver.eigenvalues()[i])<rel) rel=std::real(eigensolver.eigenvalues()[i]);
        if(std::real(eigensolver.eigenvalues()[i])>Rel) Rel=std::real(eigensolver.eigenvalues()[i]);
    }

    Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
    MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    //Shift polytope in order to contain the origin (center of the ellipsoid)
    P.shift(e);

    MT L_1 = L.inverse();
    P.linear_transformIt(L_1.transpose());

    return std::pair<NT, NT> (L_1.determinant(),rel/Rel);
}



// volume

template
<
    typename Polytope,
    typename Parameters,
    typename WalkType
>
double volume(Polytope &P,
              Parameters & var,  // constans for volume
              double error,
              unsigned int walk_length,
              WalkType walk)// = CDHRWalk<Point>())
{
    typedef typename Polytope::PolytopePoint Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> Ball;
    typedef BallIntersectPolytope<Polytope,Ball> BallPoly;
    typedef typename Parameters::RNGType RNGType;
    typedef typename Polytope::VT VT;

    //bool round = var.round;
    //bool print = var.verbose;
    //bool rand_only = var.rand_only, deltaset = false;
    unsigned int n = P.dimension();
    //unsigned int rnum = var.m;
    unsigned int rnum = std::pow(error,-2) * 400 * n * std::log(n);
    //unsigned int walk_length = var.walk_steps;
    unsigned int n_threads = 1;
    //const NT err = var.err;
    //RNGType &rng = var.rng;

    //0. Get the Chebychev ball (largest inscribed ball) with center and radius
    auto InnerBall = P.InnerBall();
    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    auto round = true;

    //1. Rounding of the polytope if round=true
    NT round_value=1;
    if(round){
#ifdef VOLESTI_DEBUG
        std::cout<<"\nRounding.."<<std::endl;
        double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
#endif
        std::pair<NT,NT> res_round = rounding_min_ellipsoid(P,InnerBall,var);
        round_value=res_round.first;
        std::pair<Point,NT> res=P.ComputeInnerBall();
        c=res.first;
        radius=res.second;
        P.comp_diam(var.diameter, radius);
        if (var.ball_walk){
            var.delta = 4.0 * radius / NT(n);
        }
#ifdef VOLESTI_DEBUG
        double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
#endif
    }

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    VT c_e = Eigen::Map<VT>(&c.get_coeffs()[0], c.dimension());
    P.shift(c_e);
    c=Point(n);

    rnum=rnum/n_threads;
    NT vol=0;

    // Perform the procedure for a number of threads and then take the average
    for(unsigned int t=0; t<n_threads; t++){
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
#ifdef VOLESTI_DEBUG
        std::cout<<"\nGenerate the first random point in P"<<std::endl;
#endif

        Point p = get_point_on_Dsphere<RNGType , Point>(n, radius);
        std::list<Point> randPoints; //ds for storing rand points

        //use a large walk length e.g. 1000
        //rand_point_generator(P, p, 1, 50*n, randPoints, var);
        PushBackWalkPolicy policy;
        walk.template apply<RNGType>(P, p, 1, 50*n, randPoints, policy);

        // 3. Sample "rnum" points from P
#ifdef VOLESTI_DEBUG
        double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
#endif
        //rand_point_generator(P, p, rnum-1, walk_len, randPoints, var);
        walk.template apply<RNGType>(P, p, rnum-1, walk_length, randPoints, policy);

#ifdef VOLESTI_DEBUG
        double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "First random points construction time = "
                  << tstop2 - tstart2 << std::endl;
#endif

        // 4.  Construct the sequence of balls
        // 4a. compute the radius of the largest ball
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
        // 4b. Number of balls
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

        // 5. Estimate Vol(P)

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

            //generate more random points in PBLarge to have "rnum" in total
            //rand_point_generator(PBLarge,p_gen,rnum-nump_PBLarge,walk_len,
            //                     randPoints,PBSmall,nump_PBSmall,var,walk);

            CountingWalkPolicy<BallPoly> counting_policy(nump_PBSmall, PBSmall);
            walk.template apply<RNGType>(PBLarge, p_gen, rnum-nump_PBLarge,
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
    vol=round_value*vol;
#ifdef VOLESTI_DEBUG
    if(print) std::cout<<"rand points = "<<rnum<<std::endl;
    if(print) std::cout<<"walk len = "<<walk_len<<std::endl;
    if(print) std::cout<<"round_value: "<<round_value<<std::endl;
    if(print) std::cout<<"volume computed: "<<vol<<std::endl;
#endif

    P.free_them_all();
    return vol;
}

#endif
