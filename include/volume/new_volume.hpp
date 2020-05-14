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

//#include "random.hpp"
//#include "random/uniform_int.hpp"
//#include "random/normal_distribution.hpp"
//#include "random/uniform_real_distribution.hpp"

#include "cartesian_geom/cartesian_kernel.h"
//#include "vars.h"
#include "new_basic_sampling_features.hpp"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "ball.h"
#include "ballintersectconvex.h"
#include "zonoIntersecthpoly.h"
#include "vpolyintersectvpoly.h"
#include "new_rounding.hpp"
//#include "samplers.h"
//#include "rounding.h"
//#include "gaussian_samplers.h"
//#include "gaussian_annealing.h"


#include "khach.h"


/////////////////// Random Walks

// ball walk with uniform target distribution
struct BallWalk
{
    BallWalk(double L)
        :   param(L, true)
    {}

    BallWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
            :   m_L(L), set_delta(set)
        {}
        double m_L;
        bool set_delta;
    };

    parameters param;

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
        typedef HPolytope<Point> Hpolytope;
        typedef Zonotope<Point> zonotope;
        typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;

        Walk (Polytope const& P, Point&, RandomNumberGenerator&) {
            _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
        }

        Walk (BallPolytope const& P, Point &, RandomNumberGenerator &) {
            _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
        }

        Walk (BallType const&, Point &, RandomNumberGenerator &) {}

        Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &) {
            _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
        }

        Walk (Polytope const& P, Point&, RandomNumberGenerator&, parameters const& params) {
            if (params.set_delta) {
                _delta = params.m_L;
            } else {
                _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
            }
        }
        Walk (BallPolytope const& P, Point &, RandomNumberGenerator &, parameters const& params) {
            if (params.set_delta) {
                _delta = params.m_L;
            } else {
                _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
            }
        }
        Walk (BallType const&, Point &, RandomNumberGenerator &, parameters &) {}

        Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &, parameters const& params) {
            if (params.set_delta) {
                _delta = params.m_L;
            } else {
                _delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());
            }
        }

        template<typename BallPolytope>
        inline void apply(BallPolytope const& P,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            //const NT delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());

            for (auto j=0u; j<walk_length; ++j)
            {
                Point y = GetPointInDsphere<Point>::apply(P.dimension(),
                                                          _delta,
                                                          rng);
                y += p;
                if (P.is_in(y) == -1) p = y;
            }
            //std::cout << "use" << parameters.m_L << std::endl;
        }

        inline void update_delta(NT delta)
        {
            _delta = delta;
        }

    private:
        double _delta;
    };

};


// random directions hit-and-run walk with uniform target distribution
struct RDHRWalk
{
    struct parameters {};
    parameters param;

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
    typedef HPolytope<Point> Hpolytope;
    typedef Zonotope<Point> zonotope;
    typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;

    Walk(Polytope const& P, Point &p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point &p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(ZonoHPoly const& P, Point &p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk (BallType const&, Point &, RandomNumberGenerator &) {}

    Walk(Polytope const& P, Point &p, RandomNumberGenerator &rng, parameters &)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point &p, RandomNumberGenerator &rng, parameters &)
    {
        initialize(P, p, rng);
    }

    Walk(ZonoHPoly const& P, Point &p, RandomNumberGenerator &rng, parameters &)
    {
        initialize(P, p, rng);
    }

    Walk (BallType const&, Point &, RandomNumberGenerator &, parameters &) {}

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
    struct parameters {};
    parameters param;

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
    typedef HPolytope<Point> Hpolytope;
    typedef Zonotope<Point> zonotope;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;
    typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;

    Walk(Polytope const& P, Point &p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point &p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk(ZonoHPoly const& P, Point &p, RandomNumberGenerator &rng)
    {
        initialize(P, p, rng);
    }

    Walk (BallType const&, Point &, RandomNumberGenerator &) {}

    Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng, parameters const& params)
    {
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point & p, RandomNumberGenerator &rng, parameters const& params)
    {
        initialize(P, p, rng);
    }

    Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &rng, parameters const& params)
    {
        initialize(P, p, rng);
    }

    Walk (BallType const&, Point &, RandomNumberGenerator &, parameters &) {}

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


struct BilliardWalkPolicy
{
    BilliardWalkPolicy(double L)
        :   _L(L)
    {}

    typedef BilliardWalk WalkType;
    double _L;
};


// random directions hit-and-run walk with uniform target distribution
struct BCDHRWalk
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

        Walk (BallType const&, Point &, RandomNumberGenerator &) {}

        template
                <
                        typename BallPolytope
                >
        inline void apply(BallPolytope const& P,
                          Point &p1,   // a point to start
                          Point &p2,
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            std::pair<NT, NT> bpair;
            for (auto j=0u; j<walk_length; ++j)
            {
                auto rand_coord_prev = _rand_coord;
                _rand_coord = rng.sample_uidist();
                NT kapa = rng.sample_urdist();
                bpair = P.line_intersect_coord(_p,
                                               _p_prev,
                                               _rand_coord,
                                               rand_coord_prev,
                                               _lamdas);
                _p_prev = _p;
                _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                                                                          * (bpair.second - bpair.first));
            }
            p1 = _p_prev;
            p2 = _p_prev;
            p1.set_coord(_rand_coord, bpair.first);
            p2.set_coord(_rand_coord, bpair.second);
            //p = _p;
        }

    private :

        template <typename GenericBody>
        inline void initialize(GenericBody const& P,
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


// random directions hit-and-run walk with uniform target distribution
struct BRDHRWalk
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
        Walk (BallType const&, Point &, RandomNumberGenerator &) {}

        template
                <
                        typename BallPolytope
                >
        inline void apply(BallPolytope const& P,
                          Point &p1,   // a point to start
                          Point &p2,
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            for (auto j=0u; j<walk_length; ++j)
            {
                Point v = GetDirection<Point>::apply(p1.dimension(), rng);
                std::pair<NT, NT> bpair = P.line_intersect(_p, v, _lamdas, _Av,
                                                           _lambda);
                _lambda = rng.sample_urdist() * (bpair.first - bpair.second)
                          + bpair.second;
                p1 = _p + bpair.first * v;
                p2 = _p + bpair.second * v;
                _p += (_lambda * v);
            }
            //p = _p;
        }

    private :

        template <typename GenericBody>
        inline void initialize(GenericBody const& P,
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


template <typename GenericPolytope>
struct compute_diameter
{
    template <typename NT>
    static NT compute(GenericPolytope) {}
};


template <typename Point>
struct compute_diameter<HPolytope<Point>>
{
template <typename NT>
static NT compute(HPolytope<Point> const& P)
{
    NT diameter = NT(4) * std::sqrt(NT(P.dimension())) * P.InnerBall().second;
    //P.set_diameter(diameter);
    return diameter;
}
};

template <typename Point>
struct compute_diameter<VPolytope<Point>>
{
template <typename NT>
static NT compute(VPolytope<Point> const& P)
{
    typedef typename VPolytope<Point>::MT MT;
    NT diameter = NT(0), diam_iter;
    MT V = P.get_mat();
    for (int i = 0; i < V.rows(); ++i) {
        for (int j = 0; j < V.rows(); ++j) {
            if (i != j) {
                diam_iter = (V.row(i) - V.row(j)).norm();
                if (diam_iter > diameter) diameter = diam_iter;
            }
        }
    }
    //P.set_diameter(diameter);
    return diameter;
}
};

template <typename Point>
struct compute_diameter<Zonotope<Point>>
{
template <typename NT>
static NT compute(Zonotope<Point> const& P)
{
    typedef typename Zonotope<Point>::MT MT;
    typedef typename Zonotope<Point>::VT VT;

    MT V = P.get_mat();
    int k = V.rows(), max_index = -1;
    MT D = V.transpose() * V;
    D = (D + D.transpose()) / 2.0;
    Eigen::SelfAdjointEigenSolver <MT> es(D);
    MT D2 = es.eigenvalues().asDiagonal(), Q = es.eigenvectors();

    NT max_eig = NT(0);
    for (int i = 0; i < P.dimension(); ++i) {
        if (es.eigenvalues()[i] > max_eig) {
            max_eig = es.eigenvalues()[i];
            max_index = i;
        }
    }

    VT max_eigvec = -1.0 * Q.col(max_index);
    VT obj_fun = max_eigvec.transpose() * V.transpose(), x0(k);

    for (int j = 0; j < k; ++j) x0(j) = (obj_fun(j) < 0.0) ? -1.0 : 1.0;

    NT diameter = NT(2) * (V.transpose() * x0).norm();
    //P.set_diameter(diameter);
    return diameter;
}
};

template <typename Point, typename RandomNumberGenerator>
struct compute_diameter<IntersectionOfVpoly<VPolytope<Point>, RandomNumberGenerator>>
{
template <typename NT>
static NT compute(IntersectionOfVpoly<VPolytope<Point>, RandomNumberGenerator> const& P)
{
    NT diameter = NT(2) * NT(P.dimension()) * P.InnerBall().second;
    //P.set_diameter(diameter);
    return diameter;
}
};

template <typename Polytope, typename Point>
struct compute_diameter<BallIntersectPolytope<Polytope, Ball<Point>>>
{
template <typename NT>
static NT compute(BallIntersectPolytope<Polytope, Ball<Point>> const& P)
{
    NT diameter = NT(2) * P.radius();
    //P.set_diameter(diameter);
    return diameter;
}
};

template <typename Point>
struct compute_diameter<ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>>>
{
template <typename NT>
static NT compute(ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>> const& P)
{
    typedef typename ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>>::VT VT;
    typedef typename ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>>::MT MT;
    typedef HPolytope<Point> Hpolytope;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    PushBackWalkPolicy push_back_policy;
    typedef typename BCDHRWalk::template Walk
            <
                    Hpolytope,
                    RandomNumberGenerator
            > BCdhrWalk;
    typedef BoundaryRandomPointGenerator<BCdhrWalk> BCdhrRandomPointGenerator;

    MT G = P.get_T().transpose();
    MT AG = P.get_mat()*G;
    int k = G.cols(), d = P.dimension();
    MT eyes1(k, 2*k);
    eyes1 << MT::Identity(k,k), NT(-1) * MT::Identity(k,k);
    MT M1(k, 4*k);
    M1 << AG.transpose(), eyes1;
    MT M = M1.transpose();
    VT b = P.get_vec();

    VT bb(4*k);
    for (int i = 0; i < 4*k; ++i) bb(i) = (i < 2*k) ? b(i) : 1.0;

    Hpolytope HP;
    HP.init(d, M, bb);

    RandomNumberGenerator rng(HP.dimension());

    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    Point q = InnerBall.first;
    BCdhrRandomPointGenerator::apply(HP, q, 4*d*d, 1,
                                     randPoints, push_back_policy, rng);
    typename std::list<Point>::iterator rpit=randPoints.begin();
    NT max_norm = NT(0), iter_norm;
    for ( ; rpit!=randPoints.end(); rpit++) {
        iter_norm = (G*(*rpit).getCoefficients()).norm();
        if (iter_norm > max_norm) max_norm = iter_norm;
    }
    NT diameter = NT(2) * max_norm;
    //P.set_diameter(diameter);
    return diameter;
}
};


// billiard walk for uniform distribution
struct BilliardWalk
{
    BilliardWalk(double L)
            :   param(L, true)
    {}

    BilliardWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

    parameters param;


template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef HPolytope<Point> Hpolytope;
    typedef Zonotope<Point> zonotope;
    typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    Walk(Polytope const& P, Point &p, RandomNumberGenerator &rng)
    {
        _L = compute_diameter<Polytope>::template compute<NT>(P);
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point &p, RandomNumberGenerator &rng)
    {
        _L = compute_diameter<BallPolytope>::template compute<NT>(P);
        initialize(P, p, rng);
    }

    Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &rng)
    {
        _L = compute_diameter<ZonoHPoly>::template compute<NT>(P);
        initialize(P, p, rng);
    }

    Walk(Polytope const& P, Point & p, RandomNumberGenerator &rng, parameters const& params)
    {
        if(params.set_L)
        {
            _L = params.m_L;
        }
        else
        {
            _L = compute_diameter<Polytope>::template compute<NT>(P);
        }
        initialize(P, p, rng);
    }

    Walk(BallPolytope const& P, Point & p, RandomNumberGenerator &rng, parameters const& params)
    {
        if(params.set_L)
        {
            _L = params.m_L;
        }
        else
        {
            _L = compute_diameter<BallPolytope>::template compute<NT>(P);
        }
        initialize(P, p, rng);
    }

    Walk(ZonoHPoly const& P, Point & p, RandomNumberGenerator &rng, parameters const& params)
    {
        if(params.set_L)
        {
            _L = params.m_L;
        }
        else
        {
            _L = compute_diameter<ZonoHPoly>::template compute<NT>(P);
        }
        initialize(P, p, rng);
    }

    Walk (BallType const&, Point &, RandomNumberGenerator &,  parameters &) {}

    Walk (BallType const&, Point &, RandomNumberGenerator &) {}


    template
    <
        typename Polytope,
        typename Point,
        typename PointList,
        typename WalkPolicy,
        typename RandomNumberGenerator,
        typename Parameters
    >
    inline void apply(GenericPolytope const& P,
                      Point &p,   // a point to start
                      unsigned int const& walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        unsigned int n = P.dimension();
        //NT diameter = P.get_diameter();
        NT T = rng.sample_urdist() * _L;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _L;
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
            if (it == 30*n) _p = p0;
        }
        p = _p;
    }

    inline void update_delta(NT L)
    {
        _L = L;
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
        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _p = p;
        _v = GetDirection<Point>::apply(n, rng);

        NT T = rng.sample_urdist() * _L;
        Point p0 = _p;
        int it = 0;

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

        while (it < 30*n)
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

    double _L;
    Point _p;
    Point _v;
    NT _lambda_prev;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
};

};


////////////////////////////// Algorithms


// ----- VOLUME ------ //

template
<
    typename WalkTypePolicy,
    typename RandomNumberGenerator,
    typename Polytope
>
double volume_sequence_of_balls(Polytope const& Pin,
                                RandomNumberGenerator &rng,
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
                                              > walk;

    typedef RandomPointGenerator<walk> RandomPointGenerator;

    auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    unsigned int rnum = std::pow(error, -2) * 400 * n * std::log(n);
    //RandomNumberGenerator rng(P.dimension());

    //Compute the Chebychev ball (largest inscribed ball) with center and radius
    auto InnerBall = P.ComputeInnerBall();
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


template
        <
                typename WalkTypePolicy,
                typename RandomNumberGenerator,
                typename Polytope
        >
double volume_sequence_of_balls(Polytope const& Pin,
                                double const& error = 1.0,
                                unsigned int const& walk_length = 1,
                                unsigned int const& n_threads = 1)
{
    RandomNumberGenerator rng(Pin.dimension());
    return volume_sequence_of_balls<WalkTypePolicy>(Pin, rng, error, walk_length, n_threads);
}


#endif
